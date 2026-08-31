from __future__ import annotations

from dataclasses import dataclass
from types import ModuleType
from typing import Any

import numpy as np

from .backend import maximum, zeros_like


@dataclass
class InfiltrationModel:
    saturated_hydraulic_conductivity: float = 0.0  # m/s
    capillary_suction: float = 0.0  # m
    porosity_deficit: float = 0.0  # [-]
    min_cumulative_depth: float = 1e-4  # m, regularization to avoid singular initial rates

    def infiltration_rate(self, h: Any, cumulative: Any, xp: ModuleType = np):
        if self.saturated_hydraulic_conductivity <= 0.0:
            return zeros_like(h, xp_module=xp)
        f = maximum(cumulative, self.min_cumulative_depth, xp_module=xp)
        green_ampt = self.saturated_hydraulic_conductivity * (1.0 + (self.capillary_suction * self.porosity_deficit) / f)
        return maximum(green_ampt, 0.0, xp_module=xp)


def apply_infiltration(
    h: Any,
    cumulative_infiltration: Any,
    dt: float,
    model: InfiltrationModel,
    xp: ModuleType = np,
):
    rate = model.infiltration_rate(h, cumulative_infiltration, xp=xp)
    infiltrated_depth = xp.minimum(rate * dt, maximum(h, 0.0, xp_module=xp))
    return h - infiltrated_depth, cumulative_infiltration + infiltrated_depth, infiltrated_depth
