from __future__ import annotations

import numpy as np


def minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))


def reconstruct_first_order(q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return q.copy(), q.copy()


def reconstruct_muscl(q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dq_minus = q - np.roll(q, 1, axis=0)
    dq_plus = np.roll(q, -1, axis=0) - q
    slope = minmod(dq_minus, dq_plus)
    q_left = q - 0.5 * slope
    q_right = q + 0.5 * slope
    return q_left, q_right
