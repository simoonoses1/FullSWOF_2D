from __future__ import annotations

from types import ModuleType
from typing import Any

import numpy as np

from .backend import maximum, sqrt


def apply_manning(hu: Any, hv: Any, h: Any, dt: float, n_manning: float, g: float, xp: ModuleType = np):
    speed = sqrt(hu * hu + hv * hv, xp_module=xp) / maximum(h, 1e-12, xp_module=xp)
    coeff = g * n_manning * n_manning * speed / maximum(h, 1e-12, xp_module=xp) ** (4.0 / 3.0)
    factor = 1.0 / (1.0 + dt * coeff)
    return hu * factor, hv * factor


def apply_darcy_weisbach(hu: Any, hv: Any, h: Any, dt: float, f_dw: float, xp: ModuleType = np):
    speed = sqrt(hu * hu + hv * hv, xp_module=xp) / maximum(h, 1e-12, xp_module=xp)
    coeff = 0.125 * f_dw * speed / maximum(h, 1e-12, xp_module=xp)
    factor = 1.0 / (1.0 + dt * coeff)
    return hu * factor, hv * factor
