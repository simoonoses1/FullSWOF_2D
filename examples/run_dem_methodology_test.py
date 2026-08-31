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
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fullswof_oil.infiltration import InfiltrationModel
from fullswof_oil.dem_workflow import ensure_interpolated_dem, interpolated_dem_path
from fullswof_oil.solver import OilSpillSolver, SolverConfig


def _optional_target_resolution(value: float | None) -> float | None:
    if value is None:
        return None
    target = float(value)
    if target <= 0.0:
        raise ValueError("--target-resolution-m must be > 0")
    return target


def _resampling_from_name(name: str) -> tuple[str, Resampling]:
    normalized = (name or "bilinear").strip().lower()
    methods = {
        "bilinear": Resampling.bilinear,
        "cubic": Resampling.cubic,
    }
    if normalized not in methods:
        raise ValueError("--dem-resampling-method must be 'bilinear' or 'cubic'")
    return normalized, methods[normalized]


def prepare_dem(
    path: Path,
    decimate: int,
    target_resolution_m: float | None = None,
    resampling_method: str = "bilinear",
) -> tuple[np.ndarray, dict, float, float]:
    target = _optional_target_resolution(target_resolution_m)
    effective_path = path
    if target is not None:
        _resampling_from_name(resampling_method)
        effective_path = ensure_interpolated_dem(path, target, resampling_method)

    with rasterio.open(effective_path) as src:
        if target is not None:
            dem = src.read(1)
            profile = src.profile
            tr = src.transform
        elif decimate <= 1:
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
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    vmin = float(np.nanmin(dem_sub))
    vmax = float(np.nanmax(dem_sub))
    if vmax <= vmin:
        vmax = vmin + 1.0
    terrain = plt.get_cmap("terrain")
    terrain_land = terrain(np.linspace(0.24, 1.0, 255))
    dark_green = np.array([[0.04, 0.22, 0.10, 1.0]])
    cmap = LinearSegmentedColormap.from_list(
        "dem_terrain_green_low",
        np.vstack([dark_green, terrain_land]),
        N=256,
    )
    return cmap, Normalize(vmin=vmin, vmax=vmax)


OIL_SPILL_CMAPS = {
    "oil_dark": [
        "#160806",
        "#2d0d07",
        "#551407",
        "#7f1d08",
        "#a52a08",
        "#c44602",
        "#d96c00",
        "#e6a12a",
    ],
    "oil_fire": [
        "#1a0b08",
        "#4a120c",
        "#8f1d14",
        "#c7361a",
        "#e85d04",
        "#f59f00",
        "#ffd166",
    ],
    "amber": [
        "#1c1006",
        "#3a1e08",
        "#63320a",
        "#8a4a0a",
        "#ad650b",
        "#cc8414",
        "#e2a83a",
    ],
}


def _build_oil_spill_cmap(palette: str = "oil_dark"):
    from matplotlib.colors import LinearSegmentedColormap

    normalized = (palette or "oil_dark").strip().lower()
    if normalized not in OIL_SPILL_CMAPS:
        valid = ", ".join(sorted(OIL_SPILL_CMAPS))
        raise ValueError(f"unknown 3D spill palette '{palette}'. Valid options: {valid}")
    return LinearSegmentedColormap.from_list(
        f"oil_spill_{normalized}",
        OIL_SPILL_CMAPS[normalized],
        N=256,
    )


def save_fig(fig, out_path: str | os.PathLike[str], dpi: int = 180) -> None:
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def crop_array(arr: np.ndarray | None, bounds):
    if arr is None:
        return None
    arr_np = np.asarray(arr)
    if bounds is None:
        return arr_np
    r0, r1, c0, c1 = bounds
    return arr_np[r0:r1, c0:c1]


def robust_vmin_vmax(arr: np.ndarray, low_pct: float = 2.0, high_pct: float = 98.0) -> tuple[float | None, float | None]:
    a = np.asarray(arr, dtype=float)
    valid = a[np.isfinite(a)]
    if valid.size == 0:
        return None, None
    pos = valid[valid > 0.0]
    sample = pos if pos.size > 0 else valid
    vmin = float(np.nanpercentile(sample, low_pct))
    vmax = float(np.nanpercentile(sample, high_pct))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        return None, None
    if vmax <= vmin:
        vmin = float(np.nanmin(sample))
        vmax = float(np.nanmax(sample))
        if vmax <= vmin:
            vmax = vmin + 1e-12
    return vmin, vmax


def _cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x, dtype=float)
    for i in range(1, len(x)):
        dt = max(float(x[i] - x[i - 1]), 0.0)
        out[i] = out[i - 1] + 0.5 * (float(y[i]) + float(y[i - 1])) * dt
    return out


def build_diag_series(outputs: list[dict], area_cell: float, q_func) -> dict[str, np.ndarray]:
    times = np.array([float(s.get("time", np.nan)) for s in outputs], dtype=float)
    q_in = np.array([max(float(q_func(float(t))), 0.0) if q_func is not None and np.isfinite(t) else 0.0 for t in times], dtype=float)
    vol_input = _cumulative_trapezoid(q_in, times)

    n = len(outputs)
    vol_surface = np.zeros(n, dtype=float)
    vol_infil = np.zeros(n, dtype=float)
    extent_m2 = np.zeros(n, dtype=float)
    centroid_row = np.full(n, np.nan, dtype=float)
    centroid_col = np.full(n, np.nan, dtype=float)

    for i, snap in enumerate(outputs):
        h = snap.get("h")
        if h is None:
            continue
        h_arr = np.asarray(h, dtype=float)
        h_pos = np.where(np.isfinite(h_arr) & (h_arr > 0.0), h_arr, 0.0)
        vol_surface[i] = float(np.sum(h_pos) * area_cell)

        wet = h_pos > 1e-9
        extent_m2[i] = float(np.count_nonzero(wet) * area_cell)

        h_sum = float(np.sum(h_pos))
        if h_sum > 0.0:
            rr, cc = np.indices(h_pos.shape)
            centroid_row[i] = float(np.sum(rr * h_pos) / h_sum)
            centroid_col[i] = float(np.sum(cc * h_pos) / h_sum)

        inf = snap.get("cumulative_infiltration")
        if inf is not None:
            inf_arr = np.asarray(inf, dtype=float)
            inf_pos = np.where(np.isfinite(inf_arr) & (inf_arr > 0.0), inf_arr, 0.0)
            vol_infil[i] = float(np.sum(inf_pos) * area_cell)

    return {
        "times": times,
        "q_in": q_in,
        "vol_input": vol_input,
        "vol_surface": vol_surface,
        "vol_infil": vol_infil,
        "extent_m2": extent_m2,
        "centroid_row": centroid_row,
        "centroid_col": centroid_col,
    }


def compute_arrival_map(outputs: list[dict], threshold: float = 1e-9) -> np.ndarray | None:
    ref = None
    for snap in outputs:
        h = snap.get("h")
        if h is not None:
            ref = np.asarray(h)
            break
    if ref is None:
        return None

    arrival = np.full(ref.shape, np.nan, dtype=float)
    for snap in outputs:
        h = snap.get("h")
        if h is None:
            continue
        h_arr = np.asarray(h, dtype=float)
        t = float(snap.get("time", np.nan))
        hit = np.isfinite(h_arr) & (h_arr > threshold) & ~np.isfinite(arrival)
        arrival[hit] = t
    return arrival


def compute_peak_maps(outputs: list[dict]) -> tuple[np.ndarray | None, np.ndarray | None]:
    ref = None
    for snap in outputs:
        h = snap.get("h")
        if h is not None:
            ref = np.asarray(h)
            break
    if ref is None:
        return None, None

    h_peak = np.full(ref.shape, -np.inf, dtype=float)
    t_peak = np.full(ref.shape, np.nan, dtype=float)
    for snap in outputs:
        h = snap.get("h")
        if h is None:
            continue
        h_arr = np.asarray(h, dtype=float)
        t = float(snap.get("time", np.nan))
        better = np.isfinite(h_arr) & (h_arr > h_peak)
        h_peak[better] = h_arr[better]
        t_peak[better] = t

    dry = ~np.isfinite(h_peak) | (h_peak <= 0.0)
    h_peak[dry] = np.nan
    t_peak[dry] = np.nan
    return h_peak, t_peak


def generate_diagnostic_report(outputs, sim, q_func, out_dir, arrival_threshold=1e-9, plot_threshold=1e-10, plot_crop_pad=8):
    if not outputs:
        return
    os.makedirs(out_dir, exist_ok=True)

    transform = None
    if hasattr(sim, "profile") and isinstance(sim.profile, dict):
        transform = sim.profile.get("transform")
    if transform is not None:
        area_cell = abs(transform.a) * abs(transform.e)
    elif hasattr(sim, "cfg"):
        area_cell = float(sim.cfg.dx * sim.cfg.dy)
    else:
        area_cell = 1.0

    series = build_diag_series(outputs, area_cell, q_func)
    crop_bounds = _crop_bounds_from_outputs(outputs, threshold=plot_threshold, pad=plot_crop_pad)
    arrival = compute_arrival_map(outputs, threshold=arrival_threshold)
    h_peak, t_peak = compute_peak_maps(outputs)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(series["times"], series["q_in"], lw=2)
    ax.set_title("Curva de caudal de entrada")
    ax.set_xlabel("Tiempo (s)")
    ax.set_ylabel("Q (m3/s)")
    ax.grid(True, alpha=0.3)
    save_fig(fig, os.path.join(out_dir, "01_curva_caudal.png"))

    inp = series["vol_input"]
    surf = series["vol_surface"]
    infil = series["vol_infil"]
    clipped = np.full_like(inp, float(getattr(sim, "V_clipped_accum", 0.0)))
    accounted = surf + infil + clipped
    diff = inp - accounted
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(series["times"], inp, label="Input acumulado", lw=2)
    ax.plot(series["times"], surf, label="Superficie", lw=1.6)
    ax.plot(series["times"], infil, label="Infiltración acumulada", lw=1.4)
    ax.plot(series["times"], clipped, label="Clipped acumulado", lw=1.4)
    ax.plot(series["times"], accounted, "--", label="Accounted", lw=1.4)
    ax.plot(series["times"], diff, ":", label="Diff", lw=1.2)
    ax.set_title("Balance de masa en el tiempo")
    ax.set_xlabel("Tiempo (s)")
    ax.set_ylabel("Volumen (m3)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    save_fig(fig, os.path.join(out_dir, "02_balance_masa.png"))

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(series["times"], series["extent_m2"], lw=2)
    ax.set_title("Área afectada (h > 1e-9 m)")
    ax.set_xlabel("Tiempo (s)")
    ax.set_ylabel("Área (m2)")
    ax.grid(True, alpha=0.3)
    save_fig(fig, os.path.join(out_dir, "03_area_afectada.png"))

    fig, ax = plt.subplots(figsize=(6, 6))
    valid = np.isfinite(series["centroid_row"]) & np.isfinite(series["centroid_col"])
    ax.plot(series["centroid_col"][valid], series["centroid_row"][valid], "-o", ms=2)
    ax.set_title("Trayectoria del centroide del derrame")
    ax.set_xlabel("Columna")
    ax.set_ylabel("Fila")
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)
    save_fig(fig, os.path.join(out_dir, "04_trayectoria_centroide.png"))

    arr_use = crop_array(arrival, crop_bounds)
    if arr_use is not None and np.isfinite(arr_use).any():
        fig, ax = plt.subplots(figsize=(8, 6))
        arr_plot = np.ma.masked_invalid(arr_use)
        im = ax.imshow(arr_plot, origin="lower", cmap="viridis")
        tmin = float(np.nanmin(arr_use))
        tmax = float(np.nanmax(arr_use))
        if tmax > tmin:
            lv = np.linspace(tmin, tmax, 12)
            cs = ax.contour(arr_plot, levels=lv, colors="white", linewidths=0.6, alpha=0.8)
            ax.clabel(cs, inline=True, fontsize=7, fmt="%.0fs")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Tiempo de llegada (s)")
        ax.set_title("Isocronas de llegada del derrame")
        ax.set_xlabel("Columna")
        ax.set_ylabel("Fila")
        save_fig(fig, os.path.join(out_dir, "05_isocronas.png"))

    hp = crop_array(h_peak, crop_bounds)
    tp = crop_array(t_peak, crop_bounds)
    if hp is not None and tp is not None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        hp_plot = np.ma.masked_where(~np.isfinite(hp) | (hp <= 0.0), hp)
        vmin, vmax = robust_vmin_vmax(hp)
        im0 = axes[0].imshow(hp_plot, origin="lower", cmap="plasma", vmin=vmin, vmax=vmax)
        axes[0].set_title("Espesor máximo (m)")
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
        tp_plot = np.ma.masked_invalid(tp)
        im1 = axes[1].imshow(tp_plot, origin="lower", cmap="magma")
        axes[1].set_title("Tiempo al pico (s)")
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
        for ax in axes:
            ax.set_xlabel("Columna")
            ax.set_ylabel("Fila")
        save_fig(fig, os.path.join(out_dir, "06_mapas_pico_tiempo.png"))

    last = outputs[-1]
    h = crop_array(last.get("h"), crop_bounds)
    zf = crop_array(last.get("zf"), crop_bounds)
    if zf is None and hasattr(sim, "z"):
        zf = crop_array(np.asarray(sim.z), crop_bounds)
    inf = crop_array(last.get("cumulative_infiltration"), crop_bounds)
    qe = crop_array(last.get("qe"), crop_bounds)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.ravel()
    datasets = [
        (h, "h final (m)", "Blues"),
        (zf, "zf final (m)", "YlOrBr"),
        (inf, "Infiltración acumulada (m)", "Greens"),
        (qe, "qe final (m/s)", "Purples"),
    ]
    for ax, (arr, title, cmap) in zip(axes, datasets):
        if arr is None:
            ax.set_axis_off()
            continue
        arr_plot = np.ma.masked_where(~np.isfinite(arr) | (arr <= 0.0), arr)
        vmin, vmax = robust_vmin_vmax(arr)
        im = ax.imshow(arr_plot, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.set_xlabel("Columna")
        ax.set_ylabel("Fila")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    save_fig(fig, os.path.join(out_dir, "07_panel_estado_final.png"))


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
    if "h" not in outputs[-1] or outputs[-1]["h"] is None:
        return
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
                  min_lift_ratio=0.006, min_lift_abs=0.25, spill_palette="oil_dark"):
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

    dem_used = dem_array[r0:r1, c0:c1] if dem_array is not None else None
    if dem_used is None:
        return
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
        h_full = snap.get("h")
        if h_full is None:
            continue
        h_crop = h_full[r0:r1, c0:c1]
        h_sub = h_crop[::step_r, ::step_c]
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
        vmin_h, vmax_h = (None, None)
    if use_log and vmin_h is not None and vmax_h is not None:
        norm_h = LogNorm(vmin=vmin_h, vmax=vmax_h)
    else:
        vmax_lin = vmax_h if vmax_h is not None and np.isfinite(vmax_h) and vmax_h > 0 else 1.0
        norm_h = Normalize(vmin=0.0, vmax=vmax_lin)

    if vmax_h is not None and np.isfinite(vmax_h) and vmax_h > 0:
        scale_units = (1000.0 if scale_mm else 1.0)
        base_thr = max(wet_threshold * scale_units, vmax_h * 1e-6)
        display_threshold_units = max(base_thr, 0.2 * vmin_h) if (vmin_h is not None and np.isfinite(vmin_h)) else base_thr
    else:
        display_threshold_units = wet_threshold * (1000.0 if scale_mm else 1.0)

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    dem_surface = ax.plot_surface(X, Y, dem_sub, cmap=dem_cmap, norm=dem_norm, linewidth=0, antialiased=False, alpha=1.0, zorder=1)
    if hasattr(dem_surface, "set_zsort"):
        dem_surface.set_zsort("min")

    spill_surface = [None]
    spill_cmap = _build_oil_spill_cmap(spill_palette)
    sm = ScalarMappable(norm=norm_h, cmap=spill_cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.08, shrink=0.8)
    cbar.set_label(f"Espesor de derrame ({'mm' if scale_mm else 'm'})")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("elevación")
    ax.view_init(elev=view_elev, azim=view_azim)
    ax.set_box_aspect((1, 1, 0.35))

    z_span = 1.0
    z_lift_min = 0.0
    if np.isfinite(dem_sub).any():
        z_dem_min = float(np.nanmin(dem_sub))
        z_dem_max = float(np.nanmax(dem_sub))
        z_span = max(z_dem_max - z_dem_min, 1e-3)
        z_lift_min = max(float(min_lift_abs), max(0.0, min_lift_ratio * z_span))
        z_extra = z_lift_min
        ax.set_zlim(z_dem_min - 0.02 * z_span, z_dem_max + max(0.1 * z_span, z_extra * 1.2))

    def _draw_spill(frame_idx):
        h_full = outputs[frame_idx].get("h")
        if h_full is None:
            return None
        h_crop = h_full[r0:r1, c0:c1]
        h_sub = h_crop[::step_r, ::step_c]
        h_plot = h_sub * (1000.0 if scale_mm else 1.0)
        wet = np.isfinite(h_plot) & (h_plot >= display_threshold_units)
        if not np.any(wet):
            return None
        z_spill = dem_sub + z_lift_min
        if hasattr(norm_h, "vmin") and norm_h.vmin is not None:
            vmin_safe = float(norm_h.vmin)
        else:
            vmin_safe = 0.0
        if hasattr(norm_h, "vmax") and norm_h.vmax is not None and np.isfinite(norm_h.vmax):
            vmax_safe = float(norm_h.vmax)
        else:
            vmax_safe = float(np.nanmax(h_plot)) if np.isfinite(h_plot).any() else max(vmin_safe, 1.0)
        if vmax_safe <= vmin_safe:
            vmax_safe = vmin_safe + 1e-12
        h_safe = np.where(h_plot > 0.0, h_plot, vmin_safe)
        h_clip = np.clip(h_safe, vmin_safe, vmax_safe)
        norm_vals = np.asarray(norm_h(h_clip), dtype=float)
        norm_vals = np.clip(norm_vals, 0.0, 1.0)
        z_spill = np.where(wet, z_spill, np.nan)
        colors = spill_cmap(norm_vals)
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
    target_resolution_m = _optional_target_resolution(args.target_resolution_m)
    dem_resampling_method = _resampling_from_name(args.dem_resampling_method)[0] if target_resolution_m is not None else None
    dem, profile, dx, dy = prepare_dem(
        Path(args.dem),
        args.decimate,
        target_resolution_m=target_resolution_m,
        resampling_method=args.dem_resampling_method,
    )
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
        friction_model=args.friction_model,
        manning_n=args.manning_n,
        darcy_f=args.darcy_f,
        fluid_density=args.fluid_density,
        dynamic_viscosity=args.dynamic_viscosity,
        yield_stress=args.yield_stress,
        h_viscous_min=args.h_viscous_min,
        evaporation_rate=args.evap_rate,
        degradation_rate=args.deg_rate,
        point_source_row=source_row,
        point_source_col=source_col,
        point_source_flow_rate=q_of_t,
        roi_buffer_m=args.roi_buffer_m,
    )
    infil = InfiltrationModel(
        saturated_hydraulic_conductivity=args.k_oil,
        capillary_suction=args.capillary_suction,
        porosity_deficit=args.porosity_deficit,
    )

    solver = OilSpillSolver(cfg, h0, z=dem, infiltration=infil)
    outputs, summary = solver.run(snapshot_interval=args.snapshot_interval)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / "final_h.npy", solver.h)
    np.save(out_dir / "final_infiltration.npy", solver.cum_infiltration)

    with open(out_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                **summary,
                "dem_stats": dem_stats,
                "dem_resampling": {
                    "decimate": args.decimate,
                    "target_resolution_m": target_resolution_m,
                    "resampling_method": dem_resampling_method,
                    "interpolated_dem_path": (
                        str(interpolated_dem_path(args.dem, target_resolution_m, args.dem_resampling_method))
                        if target_resolution_m is not None
                        else None
                    ),
                    "dx": float(dx),
                    "dy": float(dy),
                    "width": int(nx),
                    "height": int(ny),
                },
            },
            f,
            indent=2,
        )

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

    report_sim = argparse.Namespace(profile=profile, V_clipped_accum=0.0, z=dem, cfg=solver.cfg)
    generate_diagnostic_report(outputs, report_sim, q_of_t, str(out_dir))

    export_mp4(outputs, dem, profile, str(out_dir / "spill_2d.mp4"), fps=args.video_fps)
    export_mp4_3d(outputs, dem, profile, str(out_dir / "spill_3d.mp4"), fps=args.video_fps, spill_palette=args.video_3d_palette)

    print("DEM methodology run summary:")
    for k, v in summary.items():
        if isinstance(v, (int, float, np.floating)):
            print(f"  {k}: {v:.6e}")
        else:
            print(f"  {k}: {v}")
    print(f"Outputs written to: {out_dir}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run fullswof_oil methodology test over DEM")
    dem_path = r"G:\.shortcut-targets-by-id\1iuaWNgH0tpogcW8LIBph9VwoSFKajE9K\NOVA\02 - Shell\Proceso\Informes ambientales\002_MTD_ P03CdL-A\Entradas\P03CdL-A\Anal Prelim\CW502_28082024\AnalisisPreliminarLocCDL\Archivos de Ubicacion\TIF\MDT5m_Rev00.tif"
    p.add_argument("--dem", default=dem_path, help="Path to DEM GeoTIFF")
    p.add_argument("--out-dir", default="outputs_dem_test")
    p.add_argument("--decimate", type=int, default=1)
    p.add_argument("--target-resolution-m", type=float, default=None, help="Optional DEM cell size target in DEM CRS units")
    p.add_argument("--dem-resampling-method", choices=["bilinear", "cubic"], default="bilinear")
    p.add_argument("--t-end", type=float, default=1200.0)
    p.add_argument("--cfl", type=float, default=0.3)
    p.add_argument("--friction-model", choices=["none", "manning", "darcy-weisbach", "viscous"], default="manning")
    p.add_argument("--manning-n", type=float, default=0.03)
    p.add_argument("--darcy-f", type=float, default=0.02)
    p.add_argument("--fluid-density", type=float, default=900.0, help="Fluid density for viscous friction [kg/m3]")
    p.add_argument("--dynamic-viscosity", type=float, default=0.005, help="Dynamic viscosity for viscous friction [Pa*s]")
    p.add_argument("--yield-stress", type=float, default=0.0, help="Reserved yield stress parameter [Pa]")
    p.add_argument("--h-viscous-min", type=float, default=1e-4, help="Minimum wet thickness for viscous friction [m]")
    p.add_argument("--k-oil", type=float, default=5e-6)
    p.add_argument("--capillary-suction", type=float, default=0.0)
    p.add_argument("--porosity-deficit", type=float, default=0.2)
    p.add_argument("--evap-rate", type=float, default=0.0)
    p.add_argument("--deg-rate", type=float, default=0.0)
    p.add_argument("--source-row", type=int, default=None)
    p.add_argument("--source-col", type=int, default=None)
    p.add_argument("--snapshot-interval", type=float, default=10.0)
    p.add_argument("--video-fps", type=int, default=2)
    p.add_argument("--video-3d-palette", choices=sorted(OIL_SPILL_CMAPS), default="oil_dark")
    p.add_argument("--roi-buffer-m", type=float, default=500.0, help="Buffer de seguridad para el ROI dinámico (m)")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
