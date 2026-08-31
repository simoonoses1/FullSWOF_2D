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


def apply_viscous_basal(
    hu: Any,
    hv: Any,
    h: Any,
    dt: float,
    rho: float,
    mu: float,
    h_min: float = 1e-4,
    wet_tol: float | None = None,
    xp: ModuleType = np,
):
    """Semi-implicit basal viscous damping for conservative momenta.

    rho: kg/m3, mu: Pa*s = kg/(m*s), h: m, dt: s,
    lambda = 3*mu/(rho*h**2): 1/s, hu/hv: m2/s.
    """
    if rho <= 0.0:
        raise ValueError("rho must be > 0")
    if mu < 0.0:
        raise ValueError("mu must be >= 0")
    if h_min <= 0.0:
        raise ValueError("h_min must be > 0")

    threshold = h_min if wet_tol is None else wet_tol
    h_arr = xp.asarray(h, dtype=float)
    hu_new = xp.asarray(hu, dtype=float).copy()
    hv_new = xp.asarray(hv, dtype=float).copy()
    wet = xp.isfinite(h_arr) & (h_arr >= threshold)

    hu_new[~wet] = 0.0
    hv_new[~wet] = 0.0
    if mu == 0.0 or not bool(xp.any(wet)):
        return hu_new, hv_new

    h_eff = maximum(h_arr, h_min, xp_module=xp)
    damping = 1.0 + dt * 3.0 * mu / (rho * h_eff * h_eff)

    hu_new[wet] = hu_new[wet] / damping[wet]
    hv_new[wet] = hv_new[wet] / damping[wet]
    hu_new[~xp.isfinite(hu_new)] = 0.0
    hv_new[~xp.isfinite(hv_new)] = 0.0
    return hu_new, hv_new
