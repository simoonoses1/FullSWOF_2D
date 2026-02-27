from __future__ import annotations

from types import ModuleType
from typing import Any

import numpy as np

from .backend import maximum, sqrt, zeros_like

EPS = 1e-12


def _vel(hu: Any, h: Any, xp: ModuleType = np):
    out = zeros_like(hu, dtype=float, xp_module=xp)
    xp.divide(hu, h, out=out, where=h > EPS)
    return out


def physical_flux_x(h: Any, hu: Any, hv: Any, g: float, xp: ModuleType = np):
    u = _vel(hu, h, xp=xp)
    v = _vel(hv, h, xp=xp)
    return hu, hu * u + 0.5 * g * h * h, hu * v


def physical_flux_y(h: Any, hu: Any, hv: Any, g: float, xp: ModuleType = np):
    u = _vel(hu, h, xp=xp)
    v = _vel(hv, h, xp=xp)
    return hv, hv * u, hv * v + 0.5 * g * h * h


def rusanov_flux_x(
    h_l: Any,
    hu_l: Any,
    hv_l: Any,
    h_r: Any,
    hu_r: Any,
    hv_r: Any,
    g: float,
    xp: ModuleType = np,
):
    fl_h, fl_hu, fl_hv = physical_flux_x(h_l, hu_l, hv_l, g, xp=xp)
    fr_h, fr_hu, fr_hv = physical_flux_x(h_r, hu_r, hv_r, g, xp=xp)

    u_l = _vel(hu_l, h_l, xp=xp)
    u_r = _vel(hu_r, h_r, xp=xp)
    c_l = sqrt(g * maximum(h_l, 0.0, xp_module=xp), xp_module=xp)
    c_r = sqrt(g * maximum(h_r, 0.0, xp_module=xp), xp_module=xp)
    smax = maximum(xp.abs(u_l) + c_l, xp.abs(u_r) + c_r, xp_module=xp)

    f_h = 0.5 * (fl_h + fr_h) - 0.5 * smax * (h_r - h_l)
    f_hu = 0.5 * (fl_hu + fr_hu) - 0.5 * smax * (hu_r - hu_l)
    f_hv = 0.5 * (fl_hv + fr_hv) - 0.5 * smax * (hv_r - hv_l)
    return f_h, f_hu, f_hv, smax


def rusanov_flux_y(
    h_l: Any,
    hu_l: Any,
    hv_l: Any,
    h_r: Any,
    hu_r: Any,
    hv_r: Any,
    g: float,
    xp: ModuleType = np,
):
    fl_h, fl_hu, fl_hv = physical_flux_y(h_l, hu_l, hv_l, g, xp=xp)
    fr_h, fr_hu, fr_hv = physical_flux_y(h_r, hu_r, hv_r, g, xp=xp)

    v_l = _vel(hv_l, h_l, xp=xp)
    v_r = _vel(hv_r, h_r, xp=xp)
    c_l = sqrt(g * maximum(h_l, 0.0, xp_module=xp), xp_module=xp)
    c_r = sqrt(g * maximum(h_r, 0.0, xp_module=xp), xp_module=xp)
    smax = maximum(xp.abs(v_l) + c_l, xp.abs(v_r) + c_r, xp_module=xp)

    f_h = 0.5 * (fl_h + fr_h) - 0.5 * smax * (h_r - h_l)
    f_hu = 0.5 * (fl_hu + fr_hu) - 0.5 * smax * (hu_r - hu_l)
    f_hv = 0.5 * (fl_hv + fr_hv) - 0.5 * smax * (hv_r - hv_l)
    return f_h, f_hu, f_hv, smax
