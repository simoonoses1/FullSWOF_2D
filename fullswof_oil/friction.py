from __future__ import annotations

import numpy as np


def apply_manning(hu: np.ndarray, hv: np.ndarray, h: np.ndarray, dt: float, n_manning: float, g: float) -> tuple[np.ndarray, np.ndarray]:
    speed = np.sqrt(hu * hu + hv * hv) / np.maximum(h, 1e-12)
    coeff = g * n_manning * n_manning * speed / np.maximum(h, 1e-12) ** (4.0 / 3.0)
    factor = 1.0 / (1.0 + dt * coeff)
    return hu * factor, hv * factor


def apply_darcy_weisbach(hu: np.ndarray, hv: np.ndarray, h: np.ndarray, dt: float, f_dw: float) -> tuple[np.ndarray, np.ndarray]:
    speed = np.sqrt(hu * hu + hv * hv) / np.maximum(h, 1e-12)
    coeff = 0.125 * f_dw * speed / np.maximum(h, 1e-12)
    factor = 1.0 / (1.0 + dt * coeff)
    return hu * factor, hv * factor
