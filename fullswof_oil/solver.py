from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .flux import rusanov_flux_x, rusanov_flux_y
from .friction import apply_darcy_weisbach, apply_manning
from .infiltration import InfiltrationModel, apply_infiltration
from .kernels_numba import NUMBA_AVAILABLE, conservative_update_numba, max_wave_speed_numba
from .reconstruction import (
    centered_topography_source_x,
    centered_topography_source_y,
    hydrostatic_reconstruction,
)

EPS = 1e-10


def _safe_divide(num: np.ndarray, den: np.ndarray, thresh: float = EPS) -> np.ndarray:
    out = np.zeros_like(num, dtype=float)
    np.divide(num, den, out=out, where=den > thresh)
    return out


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
    point_source_row: int | None = None
    point_source_col: int | None = None
    point_source_flow_rate: Callable[[float], float] | None = None
    roi_buffer_m: float = 500.0
    execution_mode: str = "numpy"


class OilSpillSolver:
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
        if self.cfg.execution_mode not in {"numpy", "numba"}:
            raise ValueError("execution_mode must be one of {'numpy', 'numba'}")
        self.execution_mode = self.cfg.execution_mode
        if self.execution_mode == "numba" and not NUMBA_AVAILABLE:
            self.execution_mode = "numpy"

        self._cum_rain_mass = 0.0
        self._cum_evap_deg_mass = 0.0
        self._cum_point_source_mass = 0.0

        # ROI dinámico: inicializar como todo el dominio
        self.roi_slices = (slice(0, self.h.shape[0]), slice(0, self.h.shape[1]))

    def _update_roi(self, threshold: float = 1e-6):
        """Actualiza el ROI dinámico alrededor del derrame (donde h > threshold), expandido con buffer."""
        active = self.h > threshold
        if not np.any(active):
            self.roi_slices = (slice(0, self.h.shape[0]), slice(0, self.h.shape[1]))
            return
        rows, cols = np.where(active)
        min_row, max_row = rows.min(), rows.max()
        min_col, max_col = cols.min(), cols.max()
        buffer_cells_y = int(np.ceil(self.cfg.roi_buffer_m / self.cfg.dy))
        buffer_cells_x = int(np.ceil(self.cfg.roi_buffer_m / self.cfg.dx))
        roi_row_start = max(0, min_row - buffer_cells_y)
        roi_row_end   = min(self.h.shape[0], max_row + buffer_cells_y + 1)
        roi_col_start = max(0, min_col - buffer_cells_x)
        roi_col_end   = min(self.h.shape[1], max_col + buffer_cells_x + 1)
        self.roi_slices = (slice(roi_row_start, roi_row_end), slice(roi_col_start, roi_col_end))
        # ← print eliminado

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
        if self.execution_mode == "numba":
            return float(max_wave_speed_numba(self.h, self.hu, self.hv, self.cfg.g, EPS))
        u = _safe_divide(self.hu, self.h)
        v = _safe_divide(self.hv, self.h)
        c = np.sqrt(self.cfg.g * np.maximum(self.h, 0.0))
        return float(np.max(np.maximum(np.abs(u) + c, np.abs(v) + c)))

    def compute_dt(self) -> float:
        speed = self.max_wave_speed()
        if speed < 1e-14:
            return min(self.cfg.t_end - self.time, 0.1)
        dt_cfl = self.cfg.cfl * min(self.cfg.dx, self.cfg.dy) / speed
        return min(dt_cfl, self.cfg.t_end - self.time)

    def _apply_fluxes(self, dt: float) -> None:
        # Aplicar solo en el ROI
        roi = self.roi_slices
        h_roi = self.h[roi]
        hu_roi = self.hu[roi]
        hv_roi = self.hv[roi]
        z_roi = self.z[roi]

        if self.execution_mode == "numba":
            h_new, hu_new, hv_new = conservative_update_numba(
                h_roi,
                hu_roi,
                hv_roi,
                z_roi,
                dt,
                self.cfg.dx,
                self.cfg.dy,
                self.cfg.g,
                EPS,
            )
            self.h[roi] = h_new
            self.hu[roi] = hu_new
            self.hv[roi] = hv_new
            return

        hp, hup, hvp = self._apply_wall_boundaries(h_roi, hu_roi, hv_roi)
        zp = np.pad(z_roi, ((1, 1), (1, 1)), mode="edge")

        h_l, h_r = hp[1:-1, :-1], hp[1:-1, 1:]
        hu_l, hu_r = hup[1:-1, :-1], hup[1:-1, 1:]
        hv_l, hv_r = hvp[1:-1, :-1], hvp[1:-1, 1:]
        z_l, z_r = zp[1:-1, :-1], zp[1:-1, 1:]

        h_l_hr, h_r_hr = hydrostatic_reconstruction(h_l, h_r, z_l, z_r)
        u_l, u_r = _safe_divide(hu_l, h_l), _safe_divide(hu_r, h_r)
        v_l, v_r = _safe_divide(hv_l, h_l), _safe_divide(hv_r, h_r)
        fx_h, fx_hu, fx_hv, _ = rusanov_flux_x(h_l_hr, h_l_hr * u_l, h_l_hr * v_l, h_r_hr, h_r_hr * u_r, h_r_hr * v_r, self.cfg.g)
        sx = centered_topography_source_x(h_l, h_r, h_l_hr, h_r_hr, z_l, z_r, self.cfg.g)

        h_b, h_t = hp[:-1, 1:-1], hp[1:, 1:-1]
        hu_b, hu_t = hup[:-1, 1:-1], hup[1:, 1:-1]
        hv_b, hv_t = hvp[:-1, 1:-1], hvp[1:, 1:-1]
        z_b, z_t = zp[:-1, 1:-1], zp[1:, 1:-1]

        h_b_hr, h_t_hr = hydrostatic_reconstruction(h_b, h_t, z_b, z_t)
        u_b, u_t = _safe_divide(hu_b, h_b), _safe_divide(hu_t, h_t)
        v_b, v_t = _safe_divide(hv_b, h_b), _safe_divide(hv_t, h_t)
        fy_h, fy_hu, fy_hv, _ = rusanov_flux_y(h_b_hr, h_b_hr * u_b, h_b_hr * v_b, h_t_hr, h_t_hr * u_t, h_t_hr * v_t, self.cfg.g)
        sy = centered_topography_source_y(h_b, h_t, h_b_hr, h_t_hr, z_b, z_t, self.cfg.g)

        h_roi -= (dt / self.cfg.dx) * (fx_h[:, 1:] - fx_h[:, :-1]) + (dt / self.cfg.dy) * (fy_h[1:, :] - fy_h[:-1, :])
        hu_roi -= (dt / self.cfg.dx) * (fx_hu[:, 1:] - fx_hu[:, :-1]) + (dt / self.cfg.dy) * (fy_hu[1:, :] - fy_hu[:-1, :])
        hv_roi -= (dt / self.cfg.dx) * (fx_hv[:, 1:] - fx_hv[:, :-1]) + (dt / self.cfg.dy) * (fy_hv[1:, :] - fy_hv[:-1, :])

        hu_roi -= (dt / self.cfg.dx) * (sx[:, 1:] - sx[:, :-1])
        hv_roi -= (dt / self.cfg.dy) * (sy[1:, :] - sy[:-1, :])

        h_roi = np.maximum(h_roi, 0.0)
        dry = h_roi < EPS
        hu_roi[dry] = 0.0
        hv_roi[dry] = 0.0

        # Escribir de vuelta en el dominio global
        self.h[roi] = h_roi
        self.hu[roi] = hu_roi
        self.hv[roi] = hv_roi

    def _apply_sources(self, dt: float) -> None:
        roi = self.roi_slices
        cell_area = self.cfg.dx * self.cfg.dy

        # Lluvia solo en ROI
        if self.cfg.rain_rate > 0.0:
            rain_depth = self.cfg.rain_rate * dt
            self.h[roi] += rain_depth
            self._cum_rain_mass += float(rain_depth * self.cfg.nx * self.cfg.ny * cell_area)

        # Fuente puntual: si está dentro del ROI
        if self.cfg.point_source_flow_rate is not None and self.cfg.point_source_row is not None and self.cfg.point_source_col is not None:
            r = int(np.clip(self.cfg.point_source_row, 0, self.cfg.ny - 1))
            c = int(np.clip(self.cfg.point_source_col, 0, self.cfg.nx - 1))
            if (roi[0].start <= r < roi[0].stop) and (roi[1].start <= c < roi[1].stop):
                q = max(float(self.cfg.point_source_flow_rate(self.time)), 0.0)
                added_depth = q * dt / max(cell_area, 1e-12)
                self.h[r, c] += added_depth
                self._cum_point_source_mass += q * dt

        # Infiltración solo en ROI
        h_roi, cum_inf_roi, _ = apply_infiltration(self.h[roi], self.cum_infiltration[roi], dt, self.infiltration)
        self.h[roi] = h_roi
        self.cum_infiltration[roi] = cum_inf_roi

        # Fricción solo en ROI
        if self.cfg.friction_model == "manning":
            hu_roi, hv_roi = apply_manning(self.hu[roi], self.hv[roi], self.h[roi], dt, self.cfg.manning_n, self.cfg.g)
            self.hu[roi] = hu_roi
            self.hv[roi] = hv_roi
        elif self.cfg.friction_model == "darcy-weisbach":
            hu_roi, hv_roi = apply_darcy_weisbach(self.hu[roi], self.hv[roi], self.h[roi], dt, self.cfg.darcy_f)
            self.hu[roi] = hu_roi
            self.hv[roi] = hv_roi

        # Evaporación/degradación solo en ROI
        sink_rate = self.cfg.evaporation_rate + self.cfg.degradation_rate
        if sink_rate > 0.0:
            sink_depth = sink_rate * dt
            removed_depth = np.minimum(self.h[roi], sink_depth)
            self.h[roi] -= removed_depth
            self._cum_evap_deg_mass += float(np.sum(removed_depth) * cell_area)

    def step(self) -> float:
        dt = self.compute_dt()
        if dt <= 0.0:
            return 0.0

        speed_before = self.max_wave_speed()
        courant = speed_before * dt / min(self.cfg.dx, self.cfg.dy) if speed_before > 0.0 else 0.0
        self.cfl_history.append(float(courant))

        self._apply_fluxes(dt)
        self._apply_sources(dt)
        self._update_roi()
        self.time += dt
        return dt

    def run(self, snapshot_interval: float | None = None) -> tuple[list[dict], dict[str, float]]:
        outputs: list[dict] = []
        next_snap = 0.0
        step_count = 0

        def snap():
            outputs.append({
                "time": float(self.time),
                "h": self.h.copy(),
                "cumulative_infiltration": self.cum_infiltration.copy(),
            })

        initial_mass = self.total_mass
        snap()  # snapshot inicial

        while self.time < self.cfg.t_end - 1e-12:
            dt = self.step()
            if dt <= 0.0:
                break
            step_count += 1
            print(
                f"[Solver] Paso {step_count} | t={self.time:.2f} | dt={dt:.4e} | "
                f"max(h)={np.max(self.h):.4e} | "
                f"ROI filas {self.roi_slices[0].start}-{self.roi_slices[0].stop}, "
                f"cols {self.roi_slices[1].start}-{self.roi_slices[1].stop}",
                flush=True,
            )
            if snapshot_interval is not None and self.time + 1e-12 >= next_snap + snapshot_interval:
                snap()
                next_snap += snapshot_interval

        if outputs[-1]["time"] < self.time:
            snap()

        final_mass = self.total_mass
        infil_mass = float(np.sum(self.cum_infiltration) * self.cfg.dx * self.cfg.dy)
        total_inputs = initial_mass + self._cum_rain_mass + self._cum_point_source_mass
        closure = final_mass + infil_mass + self._cum_evap_deg_mass - total_inputs
        summary = {
            "initial_surface_mass": initial_mass,
            "final_surface_mass": final_mass,
            "infiltrated_mass": infil_mass,
            "evaporation_degradation_mass": self._cum_evap_deg_mass,
            "rain_input_mass": self._cum_rain_mass,
            "point_source_input_mass": self._cum_point_source_mass,
            "mass_closure_error_pct": 100.0 * closure / max(total_inputs, 1e-12),
            "max_cfl": max(self.cfl_history) if self.cfl_history else 0.0,
        }
        return outputs, summary

    @property
    def total_mass(self) -> float:
        return float(np.sum(self.h) * self.cfg.dx * self.cfg.dy)


def run_simulation(config: SolverConfig, h0: np.ndarray, z: np.ndarray | None = None, infiltration: InfiltrationModel | None = None, snapshot_interval: float | None = None):
    solver = OilSpillSolver(config, h0, z=z, infiltration=infiltration)
    outputs, summary = solver.run(snapshot_interval=snapshot_interval)
    return solver, outputs, summary
