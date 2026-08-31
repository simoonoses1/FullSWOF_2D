from __future__ import annotations

import sys
from pathlib import Path

# Allow running this script directly without installing the package.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib.pyplot as plt
import numpy as np

from fullswof_oil.infiltration import InfiltrationModel
from fullswof_oil.io import save_raster_npy
from fullswof_oil.solver import SolverConfig, run_simulation


def circular_patch(ny: int, nx: int, radius_cells: float, thickness: float) -> np.ndarray:
    y, x = np.indices((ny, nx))
    cy, cx = ny / 2.0, nx / 2.0
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    h = np.zeros((ny, nx), dtype=float)
    h[r <= radius_cells] = thickness
    return h


def main() -> None:
    nx, ny = 100, 50
    dx = dy = 1.0
    h0 = circular_patch(ny, nx, radius_cells=8.0, thickness=0.05)
    dem = np.zeros_like(h0)

    config = SolverConfig(
        nx=nx,
        ny=ny,
        dx=dx,
        dy=dy,
        t_end=1000.0,
        cfl=0.45,
        friction_model="manning",
        manning_n=0.08,
        evaporation_rate=1e-7,
        degradation_rate=1e-7,
        rain_rate=0.0,
    )
    infiltration = InfiltrationModel(
        saturated_hydraulic_conductivity=2e-6,
        capillary_suction=0.03,
        porosity_deficit=0.25,
    )
    solver, _, summary = run_simulation(config, h0, z=dem, infiltration=infiltration)

    out = REPO_ROOT / "outputs"
    out.mkdir(exist_ok=True)
    save_raster_npy(out / "final_oil_thickness.npy", solver.h)

    plt.figure(figsize=(9, 4))
    plt.imshow(solver.h, origin="lower", cmap="inferno")
    plt.colorbar(label="Oil thickness [m]")
    plt.title("Final oil thickness after 1000 s")
    plt.tight_layout()
    plt.savefig(out / "final_oil_thickness.png", dpi=150)

    print("Mass balance summary:")
    for k, v in summary.items():
        if isinstance(v, (int, float, np.floating)):
            print(f"  {k}: {v:.6e}")
        else:
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
