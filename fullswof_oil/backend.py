from __future__ import annotations

from types import ModuleType
from typing import Any

import numpy as np

try:
    import cupy as cp
except Exception:  # pragma: no cover - optional dependency
    cp = None

xp: ModuleType = np


def resolve_backend(name: str) -> ModuleType:
    backend = name.lower()
    if backend == "numpy":
        return np
    if backend == "cupy":
        if cp is None:
            raise ValueError("compute_backend='cupy' requires CuPy to be installed")
        return cp
    raise ValueError("compute_backend must be one of {'numpy', 'cupy'}")


def set_backend(name: str) -> ModuleType:
    global xp
    xp = resolve_backend(name)
    return xp


def _module(xp_module: ModuleType | None = None) -> ModuleType:
    return xp if xp_module is None else xp_module


def asarray(arr: Any, dtype: Any | None = None, xp_module: ModuleType | None = None):
    return _module(xp_module).asarray(arr, dtype=dtype)


def zeros_like(arr: Any, dtype: Any | None = None, xp_module: ModuleType | None = None):
    mod = _module(xp_module)
    if dtype is None:
        return mod.zeros_like(arr)
    return mod.zeros_like(arr, dtype=dtype)


def maximum(a: Any, b: Any, xp_module: ModuleType | None = None):
    return _module(xp_module).maximum(a, b)


def pad(arr: Any, pad_width: Any, mode: str = "constant", xp_module: ModuleType | None = None, **kwargs: Any):
    return _module(xp_module).pad(arr, pad_width, mode=mode, **kwargs)


def sqrt(arr: Any, xp_module: ModuleType | None = None):
    return _module(xp_module).sqrt(arr)


def where(condition: Any, x: Any, y: Any, xp_module: ModuleType | None = None):
    return _module(xp_module).where(condition, x, y)


def to_numpy(arr: Any) -> np.ndarray:
    if cp is not None and isinstance(arr, cp.ndarray):
        return cp.asnumpy(arr)
    return np.asarray(arr)
