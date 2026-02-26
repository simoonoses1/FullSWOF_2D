from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.enums import Resampling

from fullswof_oil.infiltration import InfiltrationModel
from fullswof_oil.solver import SolverConfig, run_simulation


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


def run(args: argparse.Namespace) -> None:
    print("Preparando DEM...")
    dem, profile, dx, dy = prepare_dem(Path(args.dem), args.decimate)
    print("DEM preparado. Validando DEM...")
    dem, dem_stats = validate_dem(dem, profile.get("nodata", None))
    print(f"DEM validado. Celdas inválidas restantes: {dem_stats['invalid_remaining']}")
    if dem_stats["invalid_remaining"] > 0:
        raise ValueError("DEM contains unresolved invalid cells")

    ny, nx = dem.shape
    print(f"Tamaño DEM: {ny} filas x {nx} columnas")
    source_row = args.source_row if args.source_row is not None else ny // 2
    source_col = args.source_col if args.source_col is not None else nx // 2

    print("Inicializando condiciones iniciales...")
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
        roi_buffer_m=args.roi_buffer_m,
    )
    infil = InfiltrationModel(
        saturated_hydraulic_conductivity=args.k_oil,
        capillary_suction=args.capillary_suction,
        porosity_deficit=args.porosity_deficit,
    )


    print("Corriendo simulación...")
    solver, summary = run_simulation(cfg, h0, z=dem, infiltration=infil)
    print("Simulación finalizada.")


    # Recortar al bounding box de celdas con flujo (h > 0)
    mask = solver.h > 0
    if np.any(mask):
        rows, cols = np.where(mask)
        min_row, max_row = rows.min(), rows.max()
        min_col, max_col = cols.min(), cols.max()
        bbox = (slice(min_row, max_row + 1), slice(min_col, max_col + 1))
        h_roi = solver.h[bbox]
        infil_roi = solver.cum_infiltration[bbox]
    else:
        h_roi = solver.h
        infil_roi = solver.cum_infiltration

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Guardando resultados en {out_dir}...")
    np.save(out_dir / "final_h.npy", h_roi)
    np.save(out_dir / "final_infiltration.npy", infil_roi)

    with open(out_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({**summary, "dem_stats": dem_stats}, f, indent=2)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    im0 = ax[0].imshow(h_roi, origin="lower", cmap="inferno")
    ax[0].set_title("Final oil thickness h [m]")
    fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)

    im1 = ax[1].imshow(infil_roi, origin="lower", cmap="Blues")
    ax[1].set_title("Cumulative infiltration [m]")
    fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)

    fig.tight_layout()
    fig.savefig(out_dir / "diagnostic_maps.png", dpi=180)
    plt.close(fig)

    print("DEM methodology run summary:")
    for k, v in summary.items():
        print(f"  {k}: {v:.6e}")
    print(f"Outputs written to: {out_dir}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run fullswof_oil methodology test over DEM")
    dem_path = r"G:\.shortcut-targets-by-id\1iuaWNgH0tpogcW8LIBph9VwoSFKajE9K\NOVA\02 - Shell\Proceso\Informes ambientales\002_MTD_ P03CdL-A\Entradas\P03CdL-A\Anal Prelim\CW502_28082024\AnalisisPreliminarLocCDL\Archivos de Ubicacion\TIF\MDT5m_Rev00.tif"
    p.add_argument("--dem", default=dem_path, help="Path to DEM GeoTIFF")
    p.add_argument("--out-dir", default="outputs_dem_test")
    p.add_argument("--decimate", type=int, default=1)
    p.add_argument("--t-end", type=float, default=150.0)
    p.add_argument("--cfl", type=float, default=0.3)
    p.add_argument("--manning-n", type=float, default=0.03)
    p.add_argument("--k-oil", type=float, default=5e-6)
    p.add_argument("--capillary-suction", type=float, default=0.0)
    p.add_argument("--porosity-deficit", type=float, default=0.2)
    p.add_argument("--evap-rate", type=float, default=0.0)
    p.add_argument("--deg-rate", type=float, default=0.0)
    p.add_argument("--source-row", type=int, default=None)
    p.add_argument("--source-col", type=int, default=None)
    p.add_argument("--roi-buffer-m", type=float, default=500.0, help="Buffer de seguridad para el ROI dinámico (m)")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
