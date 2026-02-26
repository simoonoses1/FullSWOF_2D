[1mdiff --git a/Examples/run_example.py b/Examples/run_example.py[m
[1mindex dbf2307..0dcc885 100644[m
[1m--- a/Examples/run_example.py[m
[1m+++ b/Examples/run_example.py[m
[36m@@ -1,13 +1,71 @@[m
 from __future__ import annotations[m
 [m
[31m-import runpy[m
 import sys[m
 from pathlib import Path[m
 [m
[31m-# Compatibility launcher for users running from `Examples/` (Windows/CMD, etc.)[m
[32m+[m[32m# Allow running this script directly without installing the package.[m
 REPO_ROOT = Path(__file__).resolve().parents[1][m
[31m-SCRIPT = REPO_ROOT / "examples" / "run_example.py"[m
 if str(REPO_ROOT) not in sys.path:[m
     sys.path.insert(0, str(REPO_ROOT))[m
 [m
[31m-runpy.run_path(str(SCRIPT), run_name="__main__")[m
[32m+[m[32mimport matplotlib.pyplot as plt[m
[32m+[m[32mimport numpy as np[m
[32m+[m
[32m+[m[32mfrom fullswof_oil.infiltration import InfiltrationModel[m
[32m+[m[32mfrom fullswof_oil.io import save_raster_npy[m
[32m+[m[32mfrom fullswof_oil.solver import SolverConfig, run_simulation[m
[32m+[m
[32m+[m
[32m+[m[32mdef circular_patch(ny: int, nx: int, radius_cells: float, thickness: float) -> np.ndarray:[m
[32m+[m[32m    y, x = np.indices((ny, nx))[m
[32m+[m[32m    cy, cx = ny / 2.0, nx / 2.0[m
[32m+[m[32m    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)[m
[32m+[m[32m    h = np.zeros((ny, nx), dtype=float)[m
[32m+[m[32m    h[r <= radius_cells] = thickness[m
[32m+[m[32m    return h[m
[32m+[m
[32m+[m
[32m+[m[32mdef main() -> None:[m
[32m+[m[32m    nx, ny = 100, 50[m
[32m+[m[32m    dx = dy = 1.0[m
[32m+[m[32m    h0 = circular_patch(ny, nx, radius_cells=8.0, thickness=0.05)[m
[32m+[m[32m    dem = np.zeros_like(h0)[m
[32m+[m
[32m+[m[32m    config = SolverConfig([m
[32m+[m[32m        nx=nx,[m
[32m+[m[32m        ny=ny,[m
[32m+[m[32m        dx=dx,[m
[32m+[m[32m        dy=dy,[m
[32m+[m[32m        t_end=1000.0,[m
[32m+[m[32m        cfl=0.45,[m
[32m+[m[32m        friction_model="manning",[m
[32m+[m[32m        manning_n=0.08,[m
[32m+[m[32m        evaporation_rate=1e-7,[m
[32m+[m[32m        degradation_rate=1e-7,[m
[32m+[m[32m        rain_rate=0.0,[m
[32m+[m[32m    )[m
[32m+[m[32m    infiltration = InfiltrationModel([m
[32m+[m[32m        saturated_hydraulic_conductivity=2e-6,[m
[32m+[m[32m        capillary_suction=0.03,[m
[32m+[m[32m        porosity_deficit=0.25,[m
[32m+[m[32m    )[m
[32m+[m[32m    solver, summary = run_simulation(config, h0, z=dem, infiltration=infiltration)[m
[32m+[m
[32m+[m[32m    out = REPO_ROOT / "outputs"[m
[32m+[m[32m    out.mkdir(exist_ok=True)[m
[32m+[m[32m    save_raster_npy(out / "final_oil_thickness.npy", solver.h)[m
[32m+[m
[32m+[m[32m    plt.figure(figsize=(9, 4))[m
[32m+[m[32m    plt.imshow(solver.h, origin="lower", cmap="inferno")[m
[32m+[m[32m    plt.colorbar(label="Oil thickness [m]")[m
[32m+[m[32m    plt.title("Final oil thickness after 1000 s")[m
[32m+[m[32m    plt.tight_layout()[m
[32m+[m[32m    plt.savefig(out / "final_oil_thickness.png", dpi=150)[m
[32m+[m
[32m+[m[32m    print("Mass balance summary:")[m
[32m+[m[32m    for k, v in summary.items():[m
[32m+[m[32m        print(f"  {k}: {v:.6e}")[m
[32m+[m
[32m+[m
[32m+[m[32mif __name__ == "__main__":[m
[32m+[m[32m    main()[m
