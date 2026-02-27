from __future__ import annotations

from types import ModuleType
from typing import Any

import numpy as np

from .backend import maximum


EPS = 1e-12


def minmod(a: Any, b: Any, xp: ModuleType = np):
    return 0.5 * (xp.sign(a) + xp.sign(b)) * xp.minimum(xp.abs(a), xp.abs(b))


def muscl_reconstruct(cell_values: Any, axis: int, xp: ModuleType = np):
    """Second-order MUSCL reconstruction on cell-centered values.

    Returns left/right extrapolated states in each cell along the requested axis.
    """
    d_minus = cell_values - xp.roll(cell_values, 1, axis=axis)
    d_plus = xp.roll(cell_values, -1, axis=axis) - cell_values
    slope = minmod(d_minus, d_plus, xp=xp)
    q_left = cell_values - 0.5 * slope
    q_right = cell_values + 0.5 * slope
    return q_left, q_right


def hydrostatic_reconstruction(
    h_l: Any,
    h_r: Any,
    z_l: Any,
    z_r: Any,
    xp: ModuleType = np,
):
    """Audusse hydrostatic reconstruction for well-balanced topography handling."""
    dz = z_r - z_l
    h_l_hr = maximum(0.0, h_l - maximum(0.0, dz, xp_module=xp), xp_module=xp)
    h_r_hr = maximum(0.0, h_r - maximum(0.0, -dz, xp_module=xp), xp_module=xp)
    return h_l_hr, h_r_hr


def centered_topography_source_x(
    h_l: Any,
    h_r: Any,
    h_l_hr: Any,
    h_r_hr: Any,
    z_l: Any,
    z_r: Any,
    g: float,
):
    dzc = z_r - z_l
    return 0.5 * g * (
        (h_l_hr - h_l) * (h_l_hr + h_l) + (h_r - h_r_hr) * (h_r + h_r_hr) + (h_l + h_r) * dzc
    )


def centered_topography_source_y(
    h_b: Any,
    h_t: Any,
    h_b_hr: Any,
    h_t_hr: Any,
    z_b: Any,
    z_t: Any,
    g: float,
):
    dzc = z_t - z_b
    return 0.5 * g * (
        (h_b_hr - h_b) * (h_b_hr + h_b) + (h_t - h_t_hr) * (h_t + h_t_hr) + (h_b + h_t) * dzc
    )
