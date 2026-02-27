# fullswof_oil

`fullswof_oil` is a Python 3.11 research prototype inspired by FullSWOF_2D.
It preserves a finite-volume shallow-layer structure with explicit CFL control and adapts the physics to surface oil spreading over soil with infiltration.

## Package layout

- `fullswof_oil/solver.py`: time loop, CFL timestep, mass accounting
- `fullswof_oil/flux.py`: 2D Rusanov fluxes for conservative variables
- `fullswof_oil/kernels_numba.py`: optional Numba kernels for wave speed and conservative update
- `fullswof_oil/reconstruction.py`: reconstruction helpers (first-order and MUSCL helper)
- `fullswof_oil/friction.py`: Manning and Darcy–Weisbach friction closures
- `fullswof_oil/infiltration.py`: simplified Green–Ampt-like infiltration model
- `fullswof_oil/io.py`: parameter and raster I/O
- `examples/run_example.py`: 1000 s simulation with plotting and saved output
- `tests/test_mass_conservation.py`: mass conservation and CFL checks

## Quick start

```bash
python examples/run_example.py
pytest -q tests/test_mass_conservation.py
```

The example produces:
- `outputs/final_oil_thickness.npy`
- `outputs/final_oil_thickness.png`
- terminal mass-balance summary

## Modeling notes

- State variables: `h` (oil layer thickness), `hu`, `hv`
- Numerical flux: Rusanov (local Lax–Friedrichs)
- Time integration: explicit finite-volume with CFL-limited `dt`
- Sources: friction + infiltration + optional evaporation/degradation sinks

## Execution mode selector (NumPy / Numba)

`SolverConfig` now supports `execution_mode` with values:

- `"numpy"` (default): baseline vectorized implementation.
- `"numba"`: alternative path using `@njit(parallel=True, fastmath=False)` kernels for:
  - safe wave-speed estimation,
  - Rusanov x/y fluxes,
  - conservative update + topography source terms.

If Numba is unavailable, the solver automatically falls back to `"numpy"`.

## Benchmarking wall time by mode

Run:

```bash
python examples/benchmark_execution_modes.py
```

This script reports wall-time for increasing square meshes and prints speedup (`numpy_time / numba_time`).

This code is intended for sensitivity analysis and method prototyping.

## Running without installation

If you run from the repository root:

```bash
python examples/run_example.py
```

If you are inside `Examples/` (as in Windows prompt examples):

```bash
python run_example.py
```

A compatibility launcher `Examples/run_example.py` is provided and forwards execution to `examples/run_example.py` while adding the repository root to `PYTHONPATH`.


## DEM methodology test

A DEM-based methodology test script is available:

```bash
python examples/run_dem_methodology_test.py --dem /path/to/dem.tif --out-dir outputs_dem_test
```

Key behavior:
- Loads DEM from GeoTIFF (optional decimation).
- Validates/fills nodata cells locally.
- Runs `fullswof_oil` with hydrostatic FV solver over DEM.
- Uses a point source inflow `Q(t)` (default: 0.5 m3/s up to 100 s).
- Saves final `h`, cumulative infiltration, diagnostic map PNG, `summary.json`, and video visualizations (`spill_2d.mp4` + `spill_3d.mp4`, with GIF fallback).
