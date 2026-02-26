from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.enums import Resampling

from fullswof_oil.infiltration import InfiltrationModel
from fullswof_oil.solver import OilSpillSolver, SolverConfig


def prepare_dem(path: Path, decimate: int) -> tuple[np.ndarray, dict, float, float]:
    with rasterio.open(path) as src:
        if decimate <= 1:
            dem = src.read(1)
            profile = src.profile
            tr = src.transform
        else:
            out_rows = max(1, src.height // decimate)
            out_cols = max(1, src.width // decimate)
            dem = src.read(1, out_shape=(out_rows, out_cols), resampling=Resampling.average)
            tr = src.transform * src.transform.scale(src.width / out_cols, src.height / out_rows)
            profile = src.profile.copy()
            profile.update(height=out_rows, width=out_cols, transform=tr)

    dx = abs(tr.a)
    dy = abs(tr.e)
    return dem.astype(float), profile, dx, dy


def _crop_bounds_from_outputs(outputs: list[dict], threshold: float = 1e-9, pad: int = 3):
    mask_accum = None
    for snap in outputs:
        h = snap.get("h")
        if h is None:
            continue
        wet = np.asarray(h) > threshold
        mask_accum = wet if mask_accum is None else (mask_accum | wet)
    if mask_accum is None or not np.any(mask_accum):
        return None
    rows, cols = np.where(mask_accum)
    r0 = max(int(rows.min()) - pad, 0)
    r1 = int(rows.max()) + pad + 1
    c0 = max(int(cols.min()) - pad, 0)
    c1 = int(cols.max()) + pad + 1
    return r0, r1, c0, c1


def _auto_vrange(arr: np.ndarray, percentile: float = 99.5):
    a = np.asarray(arr, dtype=float)
    pos = a[np.isfinite(a) & (a > 0)]
    if pos.size == 0:
        return 1e-6, 1.0
    vmax = float(np.nanpercentile(pos, percentile))
    vmin = float(np.nanmin(pos))
    if vmax <= vmin:
        vmax = vmin * 10.0
    return max(vmin, vmax * 1e-6), vmax


def _build_dem_cmap_norm(dem_sub: np.ndarray, sea_level: float = 0.0):
    from matplotlib.colors import Normalize

    vmin = float(np.nanmin(dem_sub))
    vmax = float(np.nanmax(dem_sub))
    if vmax <= vmin:
        vmax = vmin + 1.0
    return "terrain", Normalize(vmin=vmin, vmax=vmax)


def fill_nodata_local_average(arr: np.ndarray, invalid: np.ndarray, max_iters: int = 8) -> tuple[np.ndarray, np.ndarray]:
    work = arr.copy()
    rem = invalid.copy()
    neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    for _ in range(max_iters):
        if not rem.any():
            break
        sums = np.zeros_like(work)
        cnt = np.zeros_like(work, dtype=int)
        for dr, dc in neighbors:
            src = work[max(0, dr): work.shape[0] + min(0, dr), max(0, dc): work.shape[1] + min(0, dc)]
            src_mask = rem[max(0, dr): work.shape[0] + min(0, dr), max(0, dc): work.shape[1] + min(0, dc)]
            dst_r0, dst_r1 = max(0, -dr), work.shape[0] - max(0, dr)
            dst_c0, dst_c1 = max(0, -dc), work.shape[1] - max(0, dc)
            valid = ~src_mask
            sums[dst_r0:dst_r1, dst_c0:dst_c1] += src * valid
            cnt[dst_r0:dst_r1, dst_c0:dst_c1] += valid.astype(int)
        fillable = rem & (cnt > 0)
        if not fillable.any():
            break
        work[fillable] = sums[fillable] / cnt[fillable]
        rem[fillable] = False
    return work, rem


def validate_dem(dem: np.ndarray, nodata: float | None) -> tuple[np.ndarray, dict]:
    invalid = np.isnan(dem)
    if nodata is not None:
        invalid |= dem == nodata
    total = int(invalid.sum())
    clean, rem = fill_nodata_local_average(dem, invalid)
    stats = {
        "invalid_total": total,
        "invalid_remaining": int(rem.sum()),
        "z_min": float(np.nanmin(clean)),
        "z_max": float(np.nanmax(clean)),
    }
    return clean, stats


def q_of_t(t: float) -> float:
    return 0.5 if t <= 100.0 else 0.0


def simulate_with_snapshots(solver: OilSpillSolver, snapshot_interval: float = 10.0) -> tuple[list[dict], dict]:
    outputs: list[dict] = []
    next_snap = 0.0

    def snap():
        outputs.append(
            {
                "time": float(solver.time),
                "h": solver.h.copy(),
                "cumulative_infiltration": solver.cum_infiltration.copy(),
            }
        )

    snap()
    while solver.time < solver.cfg.t_end - 1e-12:
        if solver.step() <= 0.0:
            break
        if solver.time + 1e-12 >= next_snap + snapshot_interval:
            snap()
            next_snap += snapshot_interval
    if outputs[-1]["time"] < solver.time:
        snap()

    final_mass = solver.total_mass
    infil_mass = float(np.sum(solver.cum_infiltration) * solver.cfg.dx * solver.cfg.dy)
    total_inputs = outputs[0]["h"].sum() * solver.cfg.dx * solver.cfg.dy + solver._cum_rain_mass + solver._cum_point_source_mass
    closure = final_mass + infil_mass + solver._cum_evap_deg_mass - total_inputs
    summary = {
        "initial_surface_mass": float(outputs[0]["h"].sum() * solver.cfg.dx * solver.cfg.dy),
        "final_surface_mass": final_mass,
        "infiltrated_mass": infil_mass,
        "evaporation_degradation_mass": solver._cum_evap_deg_mass,
        "rain_input_mass": solver._cum_rain_mass,
        "point_source_input_mass": solver._cum_point_source_mass,
        "mass_closure_error_pct": 100.0 * closure / max(total_inputs, 1e-12),
        "max_cfl": max(solver.cfl_history) if solver.cfl_history else 0.0,
    }
    return outputs, summary


def export_mp4(outputs, dem_array, profile, out_path, fps=2, scale_mm=True, percentile=99.5, use_log=True, crop=True, crop_pad=2, crop_threshold=1e-9):
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from matplotlib.colors import LogNorm

    if not outputs:
        return
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    bounds = _crop_bounds_from_outputs(outputs, threshold=crop_threshold, pad=crop_pad) if crop else None
    transform = profile.get("transform", None)
    extent = None
    dem_used = dem_array
    if bounds is not None:
        r0, r1, c0, c1 = bounds
        dem_used = dem_array[r0:r1, c0:c1] if dem_array is not None else None
        if transform is not None:
            xmin = transform.c + c0 * transform.a
            xmax = transform.c + c1 * transform.a
            ymax = transform.f + r0 * transform.e
            ymin = transform.f + r1 * transform.e
            extent = (xmin, xmax, ymin, ymax)
    else:
        if transform is not None:
            xmin = transform.c
            ymax = transform.f
            xmax = xmin + transform.a * dem_array.shape[1]
            ymin = ymax + transform.e * dem_array.shape[0]
            extent = (xmin, xmax, ymin, ymax)
    h_last = outputs[-1]["h"] if bounds is None else outputs[-1]["h"][r0:r1, c0:c1]
    h_sample = h_last * (1000.0 if scale_mm else 1.0)
    vmin_h, vmax_h = _auto_vrange(h_sample, percentile=percentile)
    norm_h = LogNorm(vmin=vmin_h, vmax=vmax_h) if use_log and vmin_h is not None and vmax_h is not None else None
    fig, ax = plt.subplots(figsize=(8, 6))
    if dem_used is not None:
        ax.imshow(dem_used, origin="upper", cmap="gray", extent=extent)
    h0 = outputs[0]["h"] if bounds is None else outputs[0]["h"][r0:r1, c0:c1]
    spill_im = ax.imshow(h0 * (1000.0 if scale_mm else 1.0), origin="upper", cmap="Blues", alpha=0.7, vmin=0 if norm_h is None else None, vmax=vmax_h if norm_h is None else None, norm=norm_h, extent=extent)
    cbar = fig.colorbar(spill_im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f"Espesor ({'mm' if scale_mm else 'm'})")
    ax.set_title(f"t = {outputs[0]['time']:.1f} s")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    def init():
        spill_im.set_data(h0 * (1000.0 if scale_mm else 1.0))
        ax.set_title(f"t = {outputs[0]['time']:.1f} s")
        return [spill_im]

    def animate(i):
        h_frame_full = outputs[i]["h"]
        h_frame = h_frame_full if bounds is None else h_frame_full[r0:r1, c0:c1]
        h_frame = h_frame * (1000.0 if scale_mm else 1.0)
        spill_im.set_data(h_frame)
        ax.set_title(f"t = {outputs[i]['time']:.1f} s")
        return [spill_im]

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(outputs), interval=1000 / fps, blit=True)
    try:
        writer = animation.FFMpegWriter(fps=fps, bitrate=1800)
        anim.save(out_path, writer=writer, dpi=150)
    except Exception:
        gif_path = os.path.splitext(out_path)[0] + ".gif"
        from matplotlib.animation import PillowWriter

        anim.save(gif_path, writer=PillowWriter(fps=fps))
    plt.close(fig)


def export_mp4_3d(outputs, dem_array, profile, out_path, fps=2, scale_mm=True, percentile=99.5, use_log=True,
                  crop=True, crop_pad=2, crop_threshold=1e-9, flow_vertical_exaggeration=20.0,
                  wet_threshold=1e-9, max_grid_size=220, view_elev=50, view_azim=-125,
                  min_lift_ratio=0.006, min_lift_abs=0.25):
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from matplotlib.colors import LogNorm, Normalize
    from matplotlib.cm import ScalarMappable

    if not outputs:
        return
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    bounds = _crop_bounds_from_outputs(outputs, threshold=crop_threshold, pad=crop_pad) if crop else None
    transform = profile.get("transform", None) if profile else None

    if bounds is None:
        r0, c0 = 0, 0
        r1, c1 = dem_array.shape
    else:
        r0, r1, c0, c1 = bounds

    dem_used = dem_array[r0:r1, c0:c1]
    nr, nc = dem_used.shape
    step_r = max(1, int(math.ceil(nr / max_grid_size)))
    step_c = max(1, int(math.ceil(nc / max_grid_size)))

    row_idx = np.arange(r0, r1, step_r)
    col_idx = np.arange(c0, c1, step_c)
    dem_sub = dem_array[row_idx[:, None], col_idx[None, :]]

    if transform is not None:
        C, R = np.meshgrid(col_idx, row_idx)
        X = transform.c + C * transform.a + R * transform.b
        Y = transform.f + C * transform.d + R * transform.e
    else:
        C, R = np.meshgrid(np.arange(dem_sub.shape[1]), np.arange(dem_sub.shape[0]))
        X, Y = C, R

    dem_cmap, dem_norm = _build_dem_cmap_norm(dem_sub, sea_level=0.0)

    sampled_positive = []
    for snap in outputs:
        h_sub = snap["h"][r0:r1, c0:c1][::step_r, ::step_c]
        h_units = h_sub * (1000.0 if scale_mm else 1.0)
        pos = h_units[h_units > 0]
        if pos.size:
            sampled_positive.append(pos.ravel())
    if sampled_positive:
        h_all = np.concatenate(sampled_positive)
        vmax_h = float(np.nanpercentile(h_all, percentile))
        min_pos = float(np.nanmin(h_all))
        vmin_h = max(min_pos, vmax_h * 1e-4)
    else:
        vmin_h, vmax_h = (1e-6, 1.0)
    norm_h = LogNorm(vmin=vmin_h, vmax=vmax_h) if use_log else Normalize(vmin=0.0, vmax=vmax_h)

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    dem_surface = ax.plot_surface(X, Y, dem_sub, cmap=dem_cmap, norm=dem_norm, linewidth=0, antialiased=False, alpha=1.0, zorder=1)
    if hasattr(dem_surface, "set_zsort"):
        dem_surface.set_zsort("min")

    spill_surface = [None]
    sm = ScalarMappable(norm=norm_h, cmap="plasma")
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.08, shrink=0.8)
    cbar.set_label(f"Espesor de derrame ({'mm' if scale_mm else 'm'})")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("elevación")
    ax.view_init(elev=view_elev, azim=view_azim)
    ax.set_box_aspect((1, 1, 0.35))

    z_dem_min = float(np.nanmin(dem_sub))
    z_dem_max = float(np.nanmax(dem_sub))
    z_span = max(z_dem_max - z_dem_min, 1e-3)
    z_lift_min = max(float(min_lift_abs), max(0.0, min_lift_ratio * z_span))
    ax.set_zlim(z_dem_min - 0.02 * z_span, z_dem_max + max(0.1 * z_span, z_lift_min * 1.2))

    def _draw_spill(frame_idx):
        h_sub = outputs[frame_idx]["h"][r0:r1, c0:c1][::step_r, ::step_c]
        h_plot = h_sub * (1000.0 if scale_mm else 1.0)
        wet = np.isfinite(h_plot) & (h_plot >= wet_threshold * (1000.0 if scale_mm else 1.0))
        if not np.any(wet):
            return None
        z_spill = np.where(wet, dem_sub + z_lift_min + flow_vertical_exaggeration * (h_sub if not scale_mm else h_sub), np.nan)
        h_safe = np.where(h_plot > 0.0, h_plot, vmin_h)
        h_clip = np.clip(h_safe, vmin_h, vmax_h)
        norm_vals = np.clip(np.asarray(norm_h(h_clip), dtype=float), 0.0, 1.0)
        colors = plt.cm.plasma(norm_vals)
        colors[..., 3] = np.where(wet, 1.0, 0.0)
        spill_artist = ax.plot_surface(X, Y, z_spill, facecolors=colors, linewidth=0, edgecolor="none", antialiased=True, shade=False, zorder=10)
        if hasattr(spill_artist, "set_zsort"):
            spill_artist.set_zsort("max")
        return spill_artist

    def render_frame(i):
        if spill_surface[0] is not None:
            spill_surface[0].remove()
        spill_surface[0] = _draw_spill(i)
        ax.set_title(f"Derrame sobre DEM (3D) - t = {outputs[i]['time']:.1f} s")

    try:
        writer = animation.FFMpegWriter(fps=fps, bitrate=2200)
        with writer.saving(fig, out_path, dpi=150):
            for i in range(len(outputs)):
                render_frame(i)
                writer.grab_frame()
    except Exception:
        gif_path = os.path.splitext(out_path)[0] + ".gif"
        from matplotlib.animation import PillowWriter

        writer_gif = PillowWriter(fps=fps)
        with writer_gif.saving(fig, gif_path, dpi=150):
            for i in range(len(outputs)):
                render_frame(i)
                writer_gif.grab_frame()
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    dem, profile, dx, dy = prepare_dem(Path(args.dem), args.decimate)
    dem, dem_stats = validate_dem(dem, profile.get("nodata", None))
    if dem_stats["invalid_remaining"] > 0:
        raise ValueError("DEM contains unresolved invalid cells")

    ny, nx = dem.shape
    source_row = args.source_row if args.source_row is not None else ny // 2
    source_col = args.source_col if args.source_col is not None else nx // 2

    h0 = np.zeros_like(dem)
    cfg = SolverConfig(
        nx=nx,
        ny=ny,
        dx=dx,
        dy=dy,
        t_end=args.t_end,
        cfl=args.cfl,
        friction_model="manning",
        manning_n=args.manning_n,
        evaporation_rate=args.evap_rate,
        degradation_rate=args.deg_rate,
        point_source_row=source_row,
        point_source_col=source_col,
        point_source_flow_rate=q_of_t,
    )
    infil = InfiltrationModel(
        saturated_hydraulic_conductivity=args.k_oil,
        capillary_suction=args.capillary_suction,
        porosity_deficit=args.porosity_deficit,
    )

    solver = OilSpillSolver(cfg, h0, z=dem, infiltration=infil)
    outputs, summary = simulate_with_snapshots(solver, snapshot_interval=args.snapshot_interval)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / "final_h.npy", solver.h)
    np.save(out_dir / "final_infiltration.npy", solver.cum_infiltration)

    with open(out_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({**summary, "dem_stats": dem_stats}, f, indent=2)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    im0 = ax[0].imshow(solver.h, origin="lower", cmap="inferno")
    ax[0].set_title("Final oil thickness h [m]")
    fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)

    im1 = ax[1].imshow(solver.cum_infiltration, origin="lower", cmap="Blues")
    ax[1].set_title("Cumulative infiltration [m]")
    fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)

    fig.tight_layout()
    fig.savefig(out_dir / "diagnostic_maps.png", dpi=180)
    plt.close(fig)

    export_mp4(outputs, dem, profile, str(out_dir / "spill_2d.mp4"), fps=args.video_fps)
    export_mp4_3d(outputs, dem, profile, str(out_dir / "spill_3d.mp4"), fps=args.video_fps)

    print("DEM methodology run summary:")
    for k, v in summary.items():
        print(f"  {k}: {v:.6e}")
    print(f"Outputs written to: {out_dir}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run fullswof_oil methodology test over DEM")
    p.add_argument("--dem", required=True, help="Path to DEM GeoTIFF")
    p.add_argument("--out-dir", default="outputs_dem_test")
    p.add_argument("--decimate", type=int, default=1)
    p.add_argument("--t-end", type=float, default=2000.0)
    p.add_argument("--cfl", type=float, default=0.3)
    p.add_argument("--manning-n", type=float, default=0.03)
    p.add_argument("--k-oil", type=float, default=5e-6)
    p.add_argument("--capillary-suction", type=float, default=0.0)
    p.add_argument("--porosity-deficit", type=float, default=0.2)
    p.add_argument("--evap-rate", type=float, default=0.0)
    p.add_argument("--deg-rate", type=float, default=0.0)
    p.add_argument("--source-row", type=int, default=None)
    p.add_argument("--source-col", type=int, default=None)
    p.add_argument("--snapshot-interval", type=float, default=10.0)
    p.add_argument("--video-fps", type=int, default=2)
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
