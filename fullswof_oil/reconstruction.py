from __future__ import annotations

import numpy as np


EPS = 1e-12


def minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))


def muscl_reconstruct(cell_values: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray]:
    """Second-order MUSCL reconstruction on cell-centered values.

    Returns left/right extrapolated states in each cell along the requested axis.
    """
    d_minus = cell_values - np.roll(cell_values, 1, axis=axis)
    d_plus = np.roll(cell_values, -1, axis=axis) - cell_values
    slope = minmod(d_minus, d_plus)
    q_left = cell_values - 0.5 * slope
    q_right = cell_values + 0.5 * slope
    return q_left, q_right


def hydrostatic_reconstruction(
    h_l: np.ndarray,
    h_r: np.ndarray,
    z_l: np.ndarray,
    z_r: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Audusse hydrostatic reconstruction for well-balanced topography handling."""
    dz = z_r - z_l
    h_l_hr = np.maximum(0.0, h_l - np.maximum(0.0, dz))
    h_r_hr = np.maximum(0.0, h_r - np.maximum(0.0, -dz))
    return h_l_hr, h_r_hr


def centered_topography_source_x(
    h_l: np.ndarray,
    h_r: np.ndarray,
    h_l_hr: np.ndarray,
    h_r_hr: np.ndarray,
    z_l: np.ndarray,
    z_r: np.ndarray,
    g: float,
) -> np.ndarray:
    dzc = z_r - z_l
    return 0.5 * g * (
        (h_l_hr - h_l) * (h_l_hr + h_l) + (h_r - h_r_hr) * (h_r + h_r_hr) + (h_l + h_r) * dzc
    )


def centered_topography_source_y(
    h_b: np.ndarray,
    h_t: np.ndarray,
    h_b_hr: np.ndarray,
    h_t_hr: np.ndarray,
    z_b: np.ndarray,
    z_t: np.ndarray,
    g: float,
) -> np.ndarray:
    dzc = z_t - z_b
    return 0.5 * g * (
        (h_b_hr - h_b) * (h_b_hr + h_b) + (h_t - h_t_hr) * (h_t + h_t_hr) + (h_b + h_t) * dzc
    )
