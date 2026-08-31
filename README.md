# fullswof_oil

`fullswof_oil` is a Python 3.11 research prototype inspired by FullSWOF_2D.
It preserves a finite-volume shallow-layer structure with explicit CFL control and adapts the physics to surface oil spreading over soil with infiltration.

## Package layout

- `fullswof_oil/solver.py`: time loop, CFL timestep, mass accounting
- `fullswof_oil/flux.py`: 2D Rusanov fluxes for conservative variables
- `fullswof_oil/kernels_numba.py`: optional Numba kernels for wave speed and conservative update
- `fullswof_oil/reconstruction.py`: reconstruction helpers (first-order and MUSCL helper)
- `fullswof_oil/friction.py`: Manning, Darcy-Weisbach, and semi-implicit viscous basal friction closures
- `fullswof_oil/infiltration.py`: simplified Green-Ampt-like infiltration model
- `fullswof_oil/io.py`: parameter and raster I/O
- `Examples/run_example.py`: 1000 s simulation with plotting and saved output
- `tests/test_mass_conservation.py`: mass conservation and CFL checks

## Quick start

```bash
python Examples/run_example.py
pytest -q
```

The example produces:
- `outputs/final_oil_thickness.npy`
- `outputs/final_oil_thickness.png`
- terminal mass-balance summary

## Modeling notes

- State variables: `h` (oil layer thickness), `hu`, `hv`
- Numerical flux: Rusanov (local Lax-Friedrichs)
- Time integration: explicit finite-volume with CFL-limited `dt`
- Sources: friction + infiltration + optional evaporation/degradation sinks
- Friction models: `none`, `manning`, `darcy-weisbach`, and `viscous`
- Viscous basal friction damps only momentum with `fluid_density` [kg/m3],
  `dynamic_viscosity` [Pa*s], and `h_viscous_min` [m]; it does not remove
  oil volume directly.

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
python Examples/run_example.py
```

If you are inside `Examples/`:

```bash
python run_example.py
```

## DEM methodology test

A DEM-based methodology test script is available:

```bash
python Examples/run_dem_methodology_test.py --dem /path/to/dem.tif --out-dir outputs_dem_test
```

To use the semi-implicit viscous basal friction model:

```bash
python Examples/run_dem_methodology_test.py --dem /path/to/dem.tif --friction-model viscous --dynamic-viscosity 0.05
```

For coarse DEMs, you can run on an interpolated computational grid:

```bash
python Examples/run_dem_methodology_test.py --dem /path/to/dem_30m.tif --target-resolution-m 10 --dem-resampling-method cubic
```

Key behavior:
- Loads DEM from GeoTIFF (optional decimation).
- Validates/fills nodata cells locally.
- Runs `fullswof_oil` with hydrostatic FV solver over DEM.
- Uses a point source inflow `Q(t)` (default: 0.5 m3/s up to 100 s).
- Saves final `h`, cumulative infiltration, diagnostic map PNG, `summary.json`, and video visualizations (`spill_2d.mp4` + `spill_3d.mp4`, with GIF fallback).

## Local DEM GUI

A local web GUI for multi-site DEM runs is available:

```bash
python Examples/run_dem_gui.py
```

It opens a browser at `http://127.0.0.1:8765` and lets you select a DEM, a GPKG/SHP source layer with a `Name` field, one or more source sites, source curve parameters, fluid/model properties, and output options. Coarse DEMs can be interpolated to a target computational resolution with `bilinear` or `cubic` resampling. Each selected site is run as an independent scenario on a local DEM window around the source, controlled by `ROI buffer (m)`, under its own output folder, with a root `runs_index.json`. When enabled, each site exports a GeoTIFF of maximum oil thickness as `max_thickness_m.tif` and non-overlapping arrival-time interval polygons as `isocronas.gpkg`. Optional persisted snapshots support post-run video re-export; the GUI reports their saved size per site and can disable them to reduce disk use. If a run reports that the spill touched the crop boundary, increase the ROI buffer and rerun.
