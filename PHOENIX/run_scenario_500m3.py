from __future__ import annotations

import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fullswof_oil.dem_workflow import (  # noqa: E402
    build_cases,
    load_source_sites,
    run_site_case,
    write_runs_index,
)


SCENARIO_DIR = Path(__file__).resolve().parent
DEM_PATH = SCENARIO_DIR / "DTM CORTADO.tif"
VECTOR_PATH = SCENARIO_DIR / "Punto.gpkg"
OUTPUT_ROOT = SCENARIO_DIR / "resultado_500m3_2m"


def build_payload() -> dict:
    resolution_m = 2.0
    sites = load_source_sites(
        VECTOR_PATH,
        DEM_PATH,
        target_resolution_m=resolution_m,
        resampling_method="bilinear",
    )
    return {
        "dem_path": str(DEM_PATH),
        "vector_path": str(VECTOR_PATH),
        "out_dir": str(OUTPUT_ROOT),
        "decimate": 1,
        "target_resolution_m": resolution_m,
        "resampling_method": "bilinear",
        "sites": sites,
        "source_curve": {
            "mode": "constant",
            "volume_total_m3": 500.0,
            "duration_s": 100.0,
        },
        "model_params": {
            "t_end": 1200.0,
            "cfl": 0.3,
            "friction_model": "manning",
            "manning_n": 0.03,
            "darcy_f": 0.02,
            "fluid_density": 900.0,
            "dynamic_viscosity": 0.005,
            "yield_stress": 0.0,
            "h_viscous_min": 1e-4,
            "k_oil": 5e-6,
            "capillary_suction": 0.0,
            "porosity_deficit": 0.2,
            "evap_rate": 0.0,
            "deg_rate": 0.0,
            "snapshot_interval": 10.0,
            "roi_buffer_m": 250.0,
            "isochrone_interval_s": 60.0,
        },
        "output_options": {
            "arrays": True,
            "rasters": True,
            "maps": True,
            "diagnostic_report": True,
            "videos": False,
            "isochrones": True,
            "snapshots": True,
        },
    }


def main() -> None:
    if OUTPUT_ROOT.exists():
        raise FileExistsError(f"Output directory already exists: {OUTPUT_ROOT}")

    payload = build_payload()
    records = [run_site_case(case) for case in build_cases(payload)]
    index_path = write_runs_index(OUTPUT_ROOT, records, payload)
    print(
        json.dumps(
            {"runs_index": index_path, "records": records},
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
