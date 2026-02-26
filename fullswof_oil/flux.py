from __future__ import annotations

import numpy as np

EPS = 1e-12


def _vel(hu: np.ndarray, h: np.ndarray) -> np.ndarray:
    out = np.zeros_like(hu, dtype=float)
    np.divide(hu, h, out=out, where=h > EPS)
    return out


def physical_flux_x(h: np.ndarray, hu: np.ndarray, hv: np.ndarray, g: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u = _vel(hu, h)
    v = _vel(hv, h)
    return hu, hu * u + 0.5 * g * h * h, hu * v


def physical_flux_y(h: np.ndarray, hu: np.ndarray, hv: np.ndarray, g: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u = _vel(hu, h)
    v = _vel(hv, h)
    return hv, hv * u, hv * v + 0.5 * g * h * h


def rusanov_flux_x(
    h_l: np.ndarray,
    hu_l: np.ndarray,
    hv_l: np.ndarray,
    h_r: np.ndarray,
    hu_r: np.ndarray,
    hv_r: np.ndarray,
    g: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    fl_h, fl_hu, fl_hv = physical_flux_x(h_l, hu_l, hv_l, g)
    fr_h, fr_hu, fr_hv = physical_flux_x(h_r, hu_r, hv_r, g)

    u_l = _vel(hu_l, h_l)
    u_r = _vel(hu_r, h_r)
    c_l = np.sqrt(g * np.maximum(h_l, 0.0))
    c_r = np.sqrt(g * np.maximum(h_r, 0.0))
    smax = np.maximum(np.abs(u_l) + c_l, np.abs(u_r) + c_r)

    f_h = 0.5 * (fl_h + fr_h) - 0.5 * smax * (h_r - h_l)
    f_hu = 0.5 * (fl_hu + fr_hu) - 0.5 * smax * (hu_r - hu_l)
    f_hv = 0.5 * (fl_hv + fr_hv) - 0.5 * smax * (hv_r - hv_l)
    return f_h, f_hu, f_hv, smax


def rusanov_flux_y(
    h_l: np.ndarray,
    hu_l: np.ndarray,
    hv_l: np.ndarray,
    h_r: np.ndarray,
    hu_r: np.ndarray,
    hv_r: np.ndarray,
    g: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    fl_h, fl_hu, fl_hv = physical_flux_y(h_l, hu_l, hv_l, g)
    fr_h, fr_hu, fr_hv = physical_flux_y(h_r, hu_r, hv_r, g)

    v_l = _vel(hv_l, h_l)
    v_r = _vel(hv_r, h_r)
    c_l = np.sqrt(g * np.maximum(h_l, 0.0))
    c_r = np.sqrt(g * np.maximum(h_r, 0.0))
    smax = np.maximum(np.abs(v_l) + c_l, np.abs(v_r) + c_r)

    f_h = 0.5 * (fl_h + fr_h) - 0.5 * smax * (h_r - h_l)
    f_hu = 0.5 * (fl_hu + fr_hu) - 0.5 * smax * (hu_r - hu_l)
    f_hv = 0.5 * (fl_hv + fr_hv) - 0.5 * smax * (hv_r - hv_l)
    return f_h, f_hu, f_hv, smax
