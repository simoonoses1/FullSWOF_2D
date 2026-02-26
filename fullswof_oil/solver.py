from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .flux import rusanov_flux_x, rusanov_flux_y
from .friction import apply_darcy_weisbach, apply_manning
from .infiltration import InfiltrationModel, apply_infiltration


@dataclass
class SolverConfig:
    nx: int
    ny: int
    dx: float
    dy: float
    g: float = 9.81
    cfl: float = 0.45
    t_end: float = 100.0
    friction_model: str = "none"
    manning_n: float = 0.03
    darcy_f: float = 0.02
    evaporation_rate: float = 0.0
    degradation_rate: float = 0.0


class OilSpillSolver:
    def __init__(self, config: SolverConfig, h0: np.ndarray, z: np.ndarray | None = None, infiltration: InfiltrationModel | None = None):
        self.cfg = config
        self.h = h0.astype(float).copy()
        self.hu = np.zeros_like(self.h)
        self.hv = np.zeros_like(self.h)
        self.z = np.zeros_like(self.h) if z is None else z.astype(float).copy()
        self.infiltration = infiltration or InfiltrationModel()
        self.cum_infiltration = np.zeros_like(self.h)
        self.time = 0.0
        self.cfl_history: list[float] = []

    def max_wave_speed(self) -> float:
        u = np.where(self.h > 1e-12, self.hu / self.h, 0.0)
        v = np.where(self.h > 1e-12, self.hv / self.h, 0.0)
        c = np.sqrt(self.cfg.g * np.maximum(self.h, 0.0))
        return float(np.max(np.maximum(np.abs(u) + c, np.abs(v) + c)))

    def compute_dt(self) -> float:
        speed = self.max_wave_speed()
        if speed < 1e-14:
            return min(self.cfg.t_end - self.time, 1e-1)
        dt_cfl = self.cfg.cfl * min(self.cfg.dx, self.cfg.dy) / speed
        return min(dt_cfl, self.cfg.t_end - self.time)

    def _apply_fluxes(self, dt: float) -> None:
        hL = self.h[:, :-1]
        huL = self.hu[:, :-1]
        hvL = self.hv[:, :-1]
        hR = self.h[:, 1:]
        huR = self.hu[:, 1:]
        hvR = self.hv[:, 1:]
        fx_h, fx_hu, fx_hv = rusanov_flux_x(hL, huL, hvL, hR, huR, hvR, self.cfg.g)

        hB = self.h[:-1, :]
        huB = self.hu[:-1, :]
        hvB = self.hv[:-1, :]
        hT = self.h[1:, :]
        huT = self.hu[1:, :]
        hvT = self.hv[1:, :]
        fy_h, fy_hu, fy_hv = rusanov_flux_y(hB, huB, hvB, hT, huT, hvT, self.cfg.g)

        self.h[:, 1:-1] -= dt / self.cfg.dx * (fx_h[:, 1:] - fx_h[:, :-1])
        self.hu[:, 1:-1] -= dt / self.cfg.dx * (fx_hu[:, 1:] - fx_hu[:, :-1])
        self.hv[:, 1:-1] -= dt / self.cfg.dx * (fx_hv[:, 1:] - fx_hv[:, :-1])

        self.h[1:-1, :] -= dt / self.cfg.dy * (fy_h[1:, :] - fy_h[:-1, :])
        self.hu[1:-1, :] -= dt / self.cfg.dy * (fy_hu[1:, :] - fy_hu[:-1, :])
        self.hv[1:-1, :] -= dt / self.cfg.dy * (fy_hv[1:, :] - fy_hv[:-1, :])

        self.h = np.maximum(self.h, 0.0)
        dry = self.h < 1e-10
        self.hu[dry] = 0.0
        self.hv[dry] = 0.0

    def _apply_sources(self, dt: float) -> None:
        if self.cfg.friction_model == "manning":
            self.hu, self.hv = apply_manning(self.hu, self.hv, self.h, dt, self.cfg.manning_n, self.cfg.g)
        elif self.cfg.friction_model == "darcy-weisbach":
            self.hu, self.hv = apply_darcy_weisbach(self.hu, self.hv, self.h, dt, self.cfg.darcy_f)

        self.h, self.cum_infiltration, _ = apply_infiltration(self.h, self.cum_infiltration, dt, self.infiltration)

        sink = (self.cfg.evaporation_rate + self.cfg.degradation_rate) * dt
        if sink > 0.0:
            self.h = np.maximum(self.h - sink, 0.0)

    def step(self) -> float:
        dt = self.compute_dt()
        if dt <= 0.0:
            return 0.0
        self._apply_fluxes(dt)
        self._apply_sources(dt)
        self.time += dt
        speed = self.max_wave_speed()
        courant = speed * dt / min(self.cfg.dx, self.cfg.dy) if speed > 0.0 else 0.0
        self.cfl_history.append(float(courant))
        return dt

    def run(self) -> dict[str, float]:
        initial_mass = self.total_mass
        while self.time < self.cfg.t_end - 1e-12:
            if self.step() <= 0.0:
                break
        final_mass = self.total_mass
        infil_mass = float(np.sum(self.cum_infiltration) * self.cfg.dx * self.cfg.dy)
        return {
            "initial_surface_mass": initial_mass,
            "final_surface_mass": final_mass,
            "infiltrated_mass": infil_mass,
            "mass_error_pct": 100.0 * (final_mass + infil_mass - initial_mass) / max(initial_mass, 1e-12),
            "max_cfl": max(self.cfl_history) if self.cfl_history else 0.0,
        }

    @property
    def total_mass(self) -> float:
        return float(np.sum(self.h) * self.cfg.dx * self.cfg.dy)


def run_simulation(config: SolverConfig, h0: np.ndarray, z: np.ndarray | None = None, infiltration: InfiltrationModel | None = None):
    solver = OilSpillSolver(config, h0, z=z, infiltration=infiltration)
    summary = solver.run()
    return solver, summary
