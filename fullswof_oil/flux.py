from __future__ import annotations

import numpy as np


def physical_flux_x(h: np.ndarray, hu: np.ndarray, hv: np.ndarray, g: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u = np.where(h > 1e-12, hu / h, 0.0)
    v = np.where(h > 1e-12, hv / h, 0.0)
    return hu, hu * u + 0.5 * g * h * h, hu * v


def physical_flux_y(h: np.ndarray, hu: np.ndarray, hv: np.ndarray, g: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u = np.where(h > 1e-12, hu / h, 0.0)
    v = np.where(h > 1e-12, hv / h, 0.0)
    return hv, hv * u, hv * v + 0.5 * g * h * h


def rusanov_flux_x(
    h_l: np.ndarray,
    hu_l: np.ndarray,
    hv_l: np.ndarray,
    h_r: np.ndarray,
    hu_r: np.ndarray,
    hv_r: np.ndarray,
    g: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fl = physical_flux_x(h_l, hu_l, hv_l, g)
    fr = physical_flux_x(h_r, hu_r, hv_r, g)

    u_l = np.where(h_l > 1e-12, hu_l / h_l, 0.0)
    u_r = np.where(h_r > 1e-12, hu_r / h_r, 0.0)
    c_l = np.sqrt(g * np.maximum(h_l, 0.0))
    c_r = np.sqrt(g * np.maximum(h_r, 0.0))
    smax = np.maximum(np.abs(u_l) + c_l, np.abs(u_r) + c_r)

    ql = (h_l, hu_l, hv_l)
    qr = (h_r, hu_r, hv_r)
    return tuple(0.5 * (f_l + f_r) - 0.5 * smax * (q_r - q_l) for f_l, f_r, q_l, q_r in zip(fl, fr, ql, qr))


def rusanov_flux_y(
    h_l: np.ndarray,
    hu_l: np.ndarray,
    hv_l: np.ndarray,
    h_r: np.ndarray,
    hu_r: np.ndarray,
    hv_r: np.ndarray,
    g: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fl = physical_flux_y(h_l, hu_l, hv_l, g)
    fr = physical_flux_y(h_r, hu_r, hv_r, g)

    v_l = np.where(h_l > 1e-12, hv_l / h_l, 0.0)
    v_r = np.where(h_r > 1e-12, hv_r / h_r, 0.0)
    c_l = np.sqrt(g * np.maximum(h_l, 0.0))
    c_r = np.sqrt(g * np.maximum(h_r, 0.0))
    smax = np.maximum(np.abs(v_l) + c_l, np.abs(v_r) + c_r)

    ql = (h_l, hu_l, hv_l)
    qr = (h_r, hu_r, hv_r)
    return tuple(0.5 * (f_l + f_r) - 0.5 * smax * (q_r - q_l) for f_l, f_r, q_l, q_r in zip(fl, fr, ql, qr))
