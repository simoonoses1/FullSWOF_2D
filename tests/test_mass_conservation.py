from __future__ import annotations

import numpy as np

from fullswof_oil.infiltration import InfiltrationModel
from fullswof_oil.solver import SolverConfig, run_simulation


def test_mass_conservation_no_sinks_within_one_percent() -> None:
    ny, nx = 40, 40
    h0 = np.zeros((ny, nx), dtype=float)
    h0[15:25, 15:25] = 0.03

    cfg = SolverConfig(
        nx=nx,
        ny=ny,
        dx=1.0,
        dy=1.0,
        t_end=100.0,
        cfl=0.45,
        friction_model="none",
        evaporation_rate=0.0,
        degradation_rate=0.0,
    )
    infiltration = InfiltrationModel(0.0, 0.0, 0.0)

    _, summary = run_simulation(cfg, h0, z=np.zeros_like(h0), infiltration=infiltration)
    rel_err = abs(summary["mass_error_pct"])
    assert rel_err <= 1.0, f"mass conservation error too large: {rel_err:.3f}%"


def test_cfl_stability_respected() -> None:
    ny, nx = 20, 20
    h0 = np.zeros((ny, nx), dtype=float)
    h0[8:12, 8:12] = 0.05

    cfg = SolverConfig(nx=nx, ny=ny, dx=1.0, dy=1.0, t_end=20.0, cfl=0.4)
    solver, summary = run_simulation(cfg, h0, z=np.zeros_like(h0), infiltration=InfiltrationModel())

    assert solver.cfl_history, "CFL history should not be empty"
    assert summary["max_cfl"] <= cfg.cfl * 1.05
