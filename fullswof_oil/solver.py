from __future__ import annotations

from dataclasses import dataclass
from types import ModuleType
from typing import Any, Callable

import numpy as np

from .backend import asarray, maximum, pad, resolve_backend, set_backend, sqrt, to_numpy, zeros_like
from .flux import rusanov_flux_x, rusanov_flux_y
from .friction import apply_darcy_weisbach, apply_manning, apply_viscous_basal
from .infiltration import InfiltrationModel, apply_infiltration
from .kernels_numba import NUMBA_AVAILABLE, conservative_update_numba, max_wave_speed_numba
from .reconstruction import (
    centered_topography_source_x,
    centered_topography_source_y,
    hydrostatic_reconstruction,
)

EPS = 1e-10
VALID_FRICTION_MODELS = frozenset({"none", "manning", "darcy-weisbach", "viscous"})


def _safe_divide(num: Any, den: Any, xp: ModuleType = np, thresh: float = EPS):
    out = zeros_like(num, dtype=float, xp_module=xp)
    xp.divide(num, den, out=out, where=den > thresh)
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
    fluid_density: float = 900.0
    dynamic_viscosity: float = 0.005
    yield_stress: float = 0.0
    h_viscous_min: float = 1e-4
    evaporation_rate: float = 0.0
    degradation_rate: float = 0.0
    rain_rate: float = 0.0
    point_source_row: int | None = None
    point_source_col: int | None = None
    point_source_flow_rate: Callable[[float], float] | None = None
    roi_buffer_m: float = 500.0
    compute_backend: str = "numpy"
    execution_mode: str = "numpy"

    def __post_init__(self) -> None:
        if self.friction_model not in VALID_FRICTION_MODELS:
            valid = ", ".join(sorted(VALID_FRICTION_MODELS))
            raise ValueError(f"friction_model must be one of: {valid}")
        if self.fluid_density <= 0.0:
            raise ValueError("fluid_density must be > 0")
        if self.dynamic_viscosity < 0.0:
            raise ValueError("dynamic_viscosity must be >= 0")
        if self.yield_stress < 0.0:
            raise ValueError("yield_stress must be >= 0")
        if self.h_viscous_min <= 0.0:
            raise ValueError("h_viscous_min must be > 0")


class OilSpillSolver:
    def __init__(self, config: SolverConfig, h0: Any, z: Any | None = None, infiltration: InfiltrationModel | None = None):
        self.cfg = config
        self.xp = resolve_backend(self.cfg.compute_backend)
        set_backend(self.cfg.compute_backend)

        self.h = maximum(asarray(h0, dtype=float, xp_module=self.xp).copy(), 0.0, xp_module=self.xp)
        self.max_h = self.h.copy()
        self.hu = zeros_like(self.h, xp_module=self.xp)
        self.hv = zeros_like(self.h, xp_module=self.xp)
        self.z = zeros_like(self.h, xp_module=self.xp) if z is None else asarray(z, dtype=float, xp_module=self.xp).copy()
        self.infiltration = infiltration or InfiltrationModel()
        self.cum_infiltration = zeros_like(self.h, xp_module=self.xp)
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

        self.roi_slices = (slice(0, self.h.shape[0]), slice(0, self.h.shape[1]))

    def _update_roi(self, threshold: float = 1e-6):
        """Actualiza el ROI dinámico alrededor del derrame (donde h > threshold), expandido con buffer."""
        active = self.h > threshold
        if not bool(self.xp.any(active)):
            self.roi_slices = (slice(0, self.h.shape[0]), slice(0, self.h.shape[1]))
            return
        rows, cols = self.xp.where(active)
        min_row, max_row = int(rows.min()), int(rows.max())
        min_col, max_col = int(cols.min()), int(cols.max())
        buffer_cells_y = int(np.ceil(self.cfg.roi_buffer_m / self.cfg.dy))
        buffer_cells_x = int(np.ceil(self.cfg.roi_buffer_m / self.cfg.dx))
        roi_row_start = max(0, min_row - buffer_cells_y)
        roi_row_end = min(self.h.shape[0], max_row + buffer_cells_y + 1)
        roi_col_start = max(0, min_col - buffer_cells_x)
        roi_col_end = min(self.h.shape[1], max_col + buffer_cells_x + 1)
        self.roi_slices = (slice(roi_row_start, roi_row_end), slice(roi_col_start, roi_col_end))

    def _apply_wall_boundaries(self, h: Any, hu: Any, hv: Any):
        hp = pad(h, ((1, 1), (1, 1)), mode="edge", xp_module=self.xp)
        hup = pad(hu, ((1, 1), (1, 1)), mode="edge", xp_module=self.xp)
        hvp = pad(hv, ((1, 1), (1, 1)), mode="edge", xp_module=self.xp)
        hup[:, 0] *= -1.0
        hup[:, -1] *= -1.0
        hvp[0, :] *= -1.0
        hvp[-1, :] *= -1.0
        return hp, hup, hvp

    def max_wave_speed(self) -> float:
        if self.execution_mode == "numba":
            return float(max_wave_speed_numba(self.h, self.hu, self.hv, self.cfg.g, EPS))
        u = _safe_divide(self.hu, self.h, xp=self.xp)
        v = _safe_divide(self.hv, self.h, xp=self.xp)
        c = sqrt(self.cfg.g * maximum(self.h, 0.0, xp_module=self.xp), xp_module=self.xp)
        return float(self.xp.max(self.xp.maximum(self.xp.abs(u) + c, self.xp.abs(v) + c)))

    def compute_dt(self) -> float:
        use_roi = self.roi_slices is not None
        if use_roi:
            h_roi = self.h[self.roi_slices]
            hu_roi = self.hu[self.roi_slices]
            hv_roi = self.hv[self.roi_slices]
            use_roi = h_roi.size > 0

        if use_roi:
            u = _safe_divide(hu_roi, h_roi, xp=self.xp)
            v = _safe_divide(hv_roi, h_roi, xp=self.xp)
            c = sqrt(self.cfg.g * maximum(h_roi, 0.0, xp_module=self.xp), xp_module=self.xp)
            speed = float(self.xp.max(self.xp.maximum(self.xp.abs(u) + c, self.xp.abs(v) + c)))
        else:
            speed = self.max_wave_speed()

        if speed < 1e-14:
            return min(self.cfg.t_end - self.time, 0.1)
        dt_cfl = self.cfg.cfl * min(self.cfg.dx, self.cfg.dy) / speed
        return min(dt_cfl, self.cfg.t_end - self.time)

    def _apply_fluxes(self, dt: float) -> None:
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
        zp = pad(z_roi, ((1, 1), (1, 1)), mode="edge", xp_module=self.xp)

        h_l, h_r = hp[1:-1, :-1], hp[1:-1, 1:]
        hu_l, hu_r = hup[1:-1, :-1], hup[1:-1, 1:]
        hv_l, hv_r = hvp[1:-1, :-1], hvp[1:-1, 1:]
        z_l, z_r = zp[1:-1, :-1], zp[1:-1, 1:]

        h_l_hr, h_r_hr = hydrostatic_reconstruction(h_l, h_r, z_l, z_r, xp=self.xp)
        u_l, u_r = _safe_divide(hu_l, h_l, xp=self.xp), _safe_divide(hu_r, h_r, xp=self.xp)
        v_l, v_r = _safe_divide(hv_l, h_l, xp=self.xp), _safe_divide(hv_r, h_r, xp=self.xp)
        fx_h, fx_hu, fx_hv, _ = rusanov_flux_x(h_l_hr, h_l_hr * u_l, h_l_hr * v_l, h_r_hr, h_r_hr * u_r, h_r_hr * v_r, self.cfg.g, xp=self.xp)
        sx = centered_topography_source_x(h_l, h_r, h_l_hr, h_r_hr, z_l, z_r, self.cfg.g)

        h_b, h_t = hp[:-1, 1:-1], hp[1:, 1:-1]
        hu_b, hu_t = hup[:-1, 1:-1], hup[1:, 1:-1]
        hv_b, hv_t = hvp[:-1, 1:-1], hvp[1:, 1:-1]
        z_b, z_t = zp[:-1, 1:-1], zp[1:, 1:-1]

        h_b_hr, h_t_hr = hydrostatic_reconstruction(h_b, h_t, z_b, z_t, xp=self.xp)
        u_b, u_t = _safe_divide(hu_b, h_b, xp=self.xp), _safe_divide(hu_t, h_t, xp=self.xp)
        v_b, v_t = _safe_divide(hv_b, h_b, xp=self.xp), _safe_divide(hv_t, h_t, xp=self.xp)
        fy_h, fy_hu, fy_hv, _ = rusanov_flux_y(h_b_hr, h_b_hr * u_b, h_b_hr * v_b, h_t_hr, h_t_hr * u_t, h_t_hr * v_t, self.cfg.g, xp=self.xp)
        sy = centered_topography_source_y(h_b, h_t, h_b_hr, h_t_hr, z_b, z_t, self.cfg.g)

        h_roi -= (dt / self.cfg.dx) * (fx_h[:, 1:] - fx_h[:, :-1]) + (dt / self.cfg.dy) * (fy_h[1:, :] - fy_h[:-1, :])
        hu_roi -= (dt / self.cfg.dx) * (fx_hu[:, 1:] - fx_hu[:, :-1]) + (dt / self.cfg.dy) * (fy_hu[1:, :] - fy_hu[:-1, :])
        hv_roi -= (dt / self.cfg.dx) * (fx_hv[:, 1:] - fx_hv[:, :-1]) + (dt / self.cfg.dy) * (fy_hv[1:, :] - fy_hv[:-1, :])

        hu_roi -= (dt / self.cfg.dx) * (sx[:, 1:] - sx[:, :-1])
        hv_roi -= (dt / self.cfg.dy) * (sy[1:, :] - sy[:-1, :])

        h_roi = maximum(h_roi, 0.0, xp_module=self.xp)
        dry = h_roi < EPS
        hu_roi[dry] = 0.0
        hv_roi[dry] = 0.0

        self.h[roi] = h_roi
        self.hu[roi] = hu_roi
        self.hv[roi] = hv_roi

    def _apply_sources(self, dt: float) -> None:
        roi = self.roi_slices
        cell_area = self.cfg.dx * self.cfg.dy

        if self.cfg.rain_rate > 0.0:
            rain_depth = self.cfg.rain_rate * dt
            self.h[roi] += rain_depth
            self._cum_rain_mass += float(rain_depth * self.cfg.nx * self.cfg.ny * cell_area)

        if self.cfg.point_source_flow_rate is not None and self.cfg.point_source_row is not None and self.cfg.point_source_col is not None:
            r = int(np.clip(self.cfg.point_source_row, 0, self.cfg.ny - 1))
            c = int(np.clip(self.cfg.point_source_col, 0, self.cfg.nx - 1))
            if (roi[0].start <= r < roi[0].stop) and (roi[1].start <= c < roi[1].stop):
                q = max(float(self.cfg.point_source_flow_rate(self.time)), 0.0)
                added_depth = q * dt / max(cell_area, 1e-12)
                self.h[r, c] += added_depth
                self._cum_point_source_mass += q * dt

        h_roi, cum_inf_roi, _ = apply_infiltration(self.h[roi], self.cum_infiltration[roi], dt, self.infiltration, xp=self.xp)
        self.h[roi] = h_roi
        self.cum_infiltration[roi] = cum_inf_roi

        if self.cfg.friction_model == "manning":
            hu_roi, hv_roi = apply_manning(self.hu[roi], self.hv[roi], self.h[roi], dt, self.cfg.manning_n, self.cfg.g, xp=self.xp)
            self.hu[roi] = hu_roi
            self.hv[roi] = hv_roi
        elif self.cfg.friction_model == "darcy-weisbach":
            hu_roi, hv_roi = apply_darcy_weisbach(self.hu[roi], self.hv[roi], self.h[roi], dt, self.cfg.darcy_f, xp=self.xp)
            self.hu[roi] = hu_roi
            self.hv[roi] = hv_roi
        elif self.cfg.friction_model == "viscous":
            hu_roi, hv_roi = apply_viscous_basal(
                self.hu[roi],
                self.hv[roi],
                self.h[roi],
                dt,
                self.cfg.fluid_density,
                self.cfg.dynamic_viscosity,
                self.cfg.h_viscous_min,
                xp=self.xp,
            )
            self.hu[roi] = hu_roi
            self.hv[roi] = hv_roi

        sink_rate = self.cfg.evaporation_rate + self.cfg.degradation_rate
        if sink_rate > 0.0:
            sink_depth = sink_rate * dt
            removed_depth = self.xp.minimum(self.h[roi], sink_depth)
            self.h[roi] -= removed_depth
            self._cum_evap_deg_mass += float(self.xp.sum(removed_depth) * cell_area)

    def step(self) -> float:
        dt = self.compute_dt()
        if dt <= 0.0:
            return 0.0

        speed_before = self.max_wave_speed()
        courant = speed_before * dt / min(self.cfg.dx, self.cfg.dy) if speed_before > 0.0 else 0.0
        self.cfl_history.append(float(courant))

        self._apply_fluxes(dt)
        self._apply_sources(dt)
        self.xp.maximum(self.max_h, self.h, out=self.max_h)
        self._update_roi()
        self.time += dt
        return dt

    def run(self, snapshot_interval: float | None = None, verbose: bool = True) -> tuple[list[dict], dict[str, float | str]]:
        outputs: list[dict] = []
        next_snap = 0.0
        step_count = 0

        def snap():
            outputs.append({
                "time": float(self.time),
                "h": to_numpy(self.h.copy()),
                "cumulative_infiltration": to_numpy(self.cum_infiltration.copy()),
            })

        initial_mass = self.total_mass
        snap()

        while self.time < self.cfg.t_end - 1e-12:
            dt = self.step()
            if dt <= 0.0:
                break
            step_count += 1
            if verbose:
                print(
                    f"[Solver] Paso {step_count} | t={self.time:.2f} | dt={dt:.4e} | "
                    f"max(h)={float(self.xp.max(self.h)):.4e} | "
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
        infil_mass = float(self.xp.sum(self.cum_infiltration) * self.cfg.dx * self.cfg.dy)
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
            "max_thickness_m": float(self.xp.max(self.max_h)) if self.max_h.size else 0.0,
            "friction_model": self.cfg.friction_model,
            "fluid_density": self.cfg.fluid_density,
            "dynamic_viscosity": self.cfg.dynamic_viscosity,
            "yield_stress": self.cfg.yield_stress,
            "h_viscous_min": self.cfg.h_viscous_min,
        }
        return outputs, summary

    @property
    def total_mass(self) -> float:
        return float(self.xp.sum(self.h) * self.cfg.dx * self.cfg.dy)


def run_simulation(
    config: SolverConfig,
    h0: Any,
    z: Any | None = None,
    infiltration: InfiltrationModel | None = None,
    snapshot_interval: float | None = None,
    verbose: bool = True,
):
    solver = OilSpillSolver(config, h0, z=z, infiltration=infiltration)
    outputs, summary = solver.run(snapshot_interval=snapshot_interval, verbose=verbose)
    return solver, outputs, summary
