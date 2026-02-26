from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


try:
    import yaml
except Exception:  # pragma: no cover
    yaml = None


def load_params(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if yaml is None:
        raise RuntimeError("PyYAML is required to read params.yaml")
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def save_raster_npy(path: str | Path, arr: np.ndarray) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, arr)


def save_raster_csv(path: str | Path, arr: np.ndarray) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, arr, delimiter=",")
