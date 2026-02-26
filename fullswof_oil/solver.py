from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .flux import rusanov_flux_x, rusanov_flux_y
from .friction import apply_darcy_weisbach, apply_manning
from .infiltration import InfiltrationModel, apply_infiltration
from .reconstruction import (
    centered_topography_source_x,
    centered_topography_source_y,
    hydrostatic_reconstruction,
)

EPS = 1e-10


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
    rain_rate: float = 0.0


class OilSpillSolver:
    """Finite-volume shallow-layer solver inspired by FullSWOF_2D methodology."""

    def __init__(self, config: SolverConfig, h0: np.ndarray, z: np.ndarray | None = None, infiltration: InfiltrationModel | None = None):
        self.cfg = config
        self.h = np.maximum(h0.astype(float).copy(), 0.0)
        self.hu = np.zeros_like(self.h)
        self.hv = np.zeros_like(self.h)
        self.z = np.zeros_like(self.h) if z is None else z.astype(float).copy()
        self.infiltration = infiltration or InfiltrationModel()
        self.cum_infiltration = np.zeros_like(self.h)
        self.time = 0.0
        self.cfl_history: list[float] = []

    def _apply_wall_boundaries(self, h: np.ndarray, hu: np.ndarray, hv: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        hp = np.pad(h, ((1, 1), (1, 1)), mode="edge")
        hup = np.pad(hu, ((1, 1), (1, 1)), mode="edge")
        hvp = np.pad(hv, ((1, 1), (1, 1)), mode="edge")

        hup[:, 0] *= -1.0
        hup[:, -1] *= -1.0
        hvp[0, :] *= -1.0
        hvp[-1, :] *= -1.0
        return hp, hup, hvp

    def max_wave_speed(self) -> float:
        u = np.where(self.h > EPS, self.hu / self.h, 0.0)
        v = np.where(self.h > EPS, self.hv / self.h, 0.0)
        c = np.sqrt(self.cfg.g * np.maximum(self.h, 0.0))
        return float(np.max(np.maximum(np.abs(u) + c, np.abs(v) + c)))

    def compute_dt(self) -> float:
        speed = self.max_wave_speed()
        if speed < 1e-14:
            return min(self.cfg.t_end - self.time, 0.1)
        dt_cfl = self.cfg.cfl * min(self.cfg.dx, self.cfg.dy) / speed
        return min(dt_cfl, self.cfg.t_end - self.time)

    def _apply_fluxes(self, dt: float) -> None:
        hp, hup, hvp = self._apply_wall_boundaries(self.h, self.hu, self.hv)
        zp = np.pad(self.z, ((1, 1), (1, 1)), mode="edge")

        # Vertical interfaces (x-direction fluxes)
        h_l = hp[1:-1, :-1]
        h_r = hp[1:-1, 1:]
        hu_l = hup[1:-1, :-1]
        hu_r = hup[1:-1, 1:]
        hv_l = hvp[1:-1, :-1]
        hv_r = hvp[1:-1, 1:]
        z_l = zp[1:-1, :-1]
        z_r = zp[1:-1, 1:]

        h_l_hr, h_r_hr = hydrostatic_reconstruction(h_l, h_r, z_l, z_r)
        u_l = np.where(h_l > EPS, hu_l / h_l, 0.0)
        u_r = np.where(h_r > EPS, hu_r / h_r, 0.0)
        v_l = np.where(h_l > EPS, hv_l / h_l, 0.0)
        v_r = np.where(h_r > EPS, hv_r / h_r, 0.0)

        fx_h, fx_hu, fx_hv, _ = rusanov_flux_x(h_l_hr, h_l_hr * u_l, h_l_hr * v_l, h_r_hr, h_r_hr * u_r, h_r_hr * v_r, self.cfg.g)
        sx = centered_topography_source_x(h_l, h_r, h_l_hr, h_r_hr, z_l, z_r, self.cfg.g)

        # Horizontal interfaces (y-direction fluxes)
        h_b = hp[:-1, 1:-1]
        h_t = hp[1:, 1:-1]
        hu_b = hup[:-1, 1:-1]
        hu_t = hup[1:, 1:-1]
        hv_b = hvp[:-1, 1:-1]
        hv_t = hvp[1:, 1:-1]
        z_b = zp[:-1, 1:-1]
        z_t = zp[1:, 1:-1]

        h_b_hr, h_t_hr = hydrostatic_reconstruction(h_b, h_t, z_b, z_t)
        u_b = np.where(h_b > EPS, hu_b / h_b, 0.0)
        u_t = np.where(h_t > EPS, hu_t / h_t, 0.0)
        v_b = np.where(h_b > EPS, hv_b / h_b, 0.0)
        v_t = np.where(h_t > EPS, hv_t / h_t, 0.0)

        fy_h, fy_hu, fy_hv, _ = rusanov_flux_y(h_b_hr, h_b_hr * u_b, h_b_hr * v_b, h_t_hr, h_t_hr * u_t, h_t_hr * v_t, self.cfg.g)
        sy = centered_topography_source_y(h_b, h_t, h_b_hr, h_t_hr, z_b, z_t, self.cfg.g)

        # Update interior cells
        self.h -= (dt / self.cfg.dx) * (fx_h[:, 1:] - fx_h[:, :-1]) + (dt / self.cfg.dy) * (fy_h[1:, :] - fy_h[:-1, :])
        self.hu -= (dt / self.cfg.dx) * (fx_hu[:, 1:] - fx_hu[:, :-1]) + (dt / self.cfg.dy) * (fy_hu[1:, :] - fy_hu[:-1, :])
        self.hv -= (dt / self.cfg.dx) * (fx_hv[:, 1:] - fx_hv[:, :-1]) + (dt / self.cfg.dy) * (fy_hv[1:, :] - fy_hv[:-1, :])

        # Bed-slope correction (well-balanced form)
        self.hu -= (dt / self.cfg.dx) * (sx[:, 1:] - sx[:, :-1])
        self.hv -= (dt / self.cfg.dy) * (sy[1:, :] - sy[:-1, :])

        self.h = np.maximum(self.h, 0.0)
        dry = self.h < EPS
        self.hu[dry] = 0.0
        self.hv[dry] = 0.0

    def _apply_sources(self, dt: float) -> None:
        # Rainfall input (as thickness source)
        if self.cfg.rain_rate > 0.0:
            self.h += self.cfg.rain_rate * dt

        self.h, self.cum_infiltration, _ = apply_infiltration(self.h, self.cum_infiltration, dt, self.infiltration)

        if self.cfg.friction_model == "manning":
            self.hu, self.hv = apply_manning(self.hu, self.hv, self.h, dt, self.cfg.manning_n, self.cfg.g)
        elif self.cfg.friction_model == "darcy-weisbach":
            self.hu, self.hv = apply_darcy_weisbach(self.hu, self.hv, self.h, dt, self.cfg.darcy_f)

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
        expected_sink = (self.cfg.evaporation_rate + self.cfg.degradation_rate) * self.cfg.nx * self.cfg.ny * self.cfg.dx * self.cfg.dy * self.time
        closure = final_mass + infil_mass + expected_sink - initial_mass
        return {
            "initial_surface_mass": initial_mass,
            "final_surface_mass": final_mass,
            "infiltrated_mass": infil_mass,
            "expected_uniform_sink_mass": expected_sink,
            "mass_closure_error_pct": 100.0 * closure / max(initial_mass, 1e-12),
            "max_cfl": max(self.cfl_history) if self.cfl_history else 0.0,
        }

    @property
    def total_mass(self) -> float:
        return float(np.sum(self.h) * self.cfg.dx * self.cfg.dy)


def run_simulation(config: SolverConfig, h0: np.ndarray, z: np.ndarray | None = None, infiltration: InfiltrationModel | None = None):
    solver = OilSpillSolver(config, h0, z=z, infiltration=infiltration)
    summary = solver.run()
    return solver, summary
