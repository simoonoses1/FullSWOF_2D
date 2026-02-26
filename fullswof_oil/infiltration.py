from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class InfiltrationModel:
    saturated_hydraulic_conductivity: float = 0.0  # m/s
    capillary_suction: float = 0.0  # m
    porosity_deficit: float = 0.0  # [-]
    min_cumulative_depth: float = 1e-4  # m, regularization to avoid singular initial rates

    def infiltration_rate(self, h: np.ndarray, cumulative: np.ndarray) -> np.ndarray:
        if self.saturated_hydraulic_conductivity <= 0.0:
            return np.zeros_like(h)
        f = np.maximum(cumulative, self.min_cumulative_depth)
        green_ampt = self.saturated_hydraulic_conductivity * (1.0 + (self.capillary_suction * self.porosity_deficit) / f)
        return np.maximum(green_ampt, 0.0)


def apply_infiltration(
    h: np.ndarray,
    cumulative_infiltration: np.ndarray,
    dt: float,
    model: InfiltrationModel,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rate = model.infiltration_rate(h, cumulative_infiltration)
    infiltrated_depth = np.minimum(rate * dt, np.maximum(h, 0.0))
    return h - infiltrated_depth, cumulative_infiltration + infiltrated_depth, infiltrated_depth
