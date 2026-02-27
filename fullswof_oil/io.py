from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .backend import to_numpy


def load_params(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def save_raster_npy(path: str | Path, arr) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, to_numpy(arr))


def save_raster_csv(path: str | Path, arr) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, to_numpy(arr), delimiter=",")
