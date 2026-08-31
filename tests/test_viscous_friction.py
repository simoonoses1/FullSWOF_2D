from __future__ import annotations

import numpy as np
import pytest

from fullswof_oil.friction import apply_viscous_basal
from fullswof_oil.solver import SolverConfig


def test_viscous_zero_viscosity_keeps_wet_momentum() -> None:
    h = np.array([[0.02, 0.05], [0.03, 0.01]])
    hu = np.array([[1.0, -2.0], [0.5, -0.25]])
    hv = np.array([[-0.4, 0.8], [1.2, -1.5]])

    hu_new, hv_new = apply_viscous_basal(hu, hv, h, dt=2.0, rho=900.0, mu=0.0, h_min=1e-4)

    np.testing.assert_allclose(hu_new, hu)
    np.testing.assert_allclose(hv_new, hv)


def test_viscous_positive_viscosity_dissipates_momentum() -> None:
    h = np.array([[0.02, 0.05], [0.03, 0.01]])
    hu = np.array([[1.0, -2.0], [0.5, -0.25]])
    hv = np.array([[-0.4, 0.8], [1.2, -1.5]])

    hu_new, hv_new = apply_viscous_basal(hu, hv, h, dt=2.0, rho=900.0, mu=0.02, h_min=1e-4)

    assert np.all(np.abs(hu_new) <= np.abs(hu))
    assert np.all(np.abs(hv_new) <= np.abs(hv))
    assert np.any(np.abs(hu_new) < np.abs(hu))
    assert np.any(np.abs(hv_new) < np.abs(hv))


def test_viscous_friction_preserves_momentum_sign() -> None:
    h = np.array([[0.02, 0.02], [0.02, 0.02]])
    hu = np.array([[1.0, -2.0], [0.5, -0.25]])
    hv = np.array([[-0.4, 0.8], [1.2, -1.5]])

    hu_new, hv_new = apply_viscous_basal(hu, hv, h, dt=1.0, rho=900.0, mu=0.1, h_min=1e-4)

    assert np.array_equal(np.sign(hu_new), np.sign(hu))
    assert np.array_equal(np.sign(hv_new), np.sign(hv))


def test_viscous_friction_zeros_dry_cells() -> None:
    h = np.array([[0.0, 5e-5], [1e-4, 2e-4]])
    hu = np.ones_like(h)
    hv = -np.ones_like(h)

    hu_new, hv_new = apply_viscous_basal(hu, hv, h, dt=1.0, rho=900.0, mu=0.02, h_min=1e-4)

    dry = h < 1e-4
    assert np.all(hu_new[dry] == 0.0)
    assert np.all(hv_new[dry] == 0.0)
    assert np.all(hu_new[~dry] > 0.0)
    assert np.all(hv_new[~dry] < 0.0)


def test_viscous_friction_matches_local_analytic_update() -> None:
    h = np.array([[0.03]])
    hu = np.array([[1.4]])
    hv = np.array([[-0.7]])
    dt = 3.0
    rho = 850.0
    mu = 0.015

    hu_new, hv_new = apply_viscous_basal(hu, hv, h, dt=dt, rho=rho, mu=mu, h_min=1e-4)

    lam = 3.0 * mu / (rho * h[0, 0] ** 2)
    expected_factor = 1.0 / (1.0 + dt * lam)
    assert np.isclose(hu_new[0, 0], hu[0, 0] * expected_factor)
    assert np.isclose(hv_new[0, 0], hv[0, 0] * expected_factor)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"friction_model": "bogus"},
        {"fluid_density": 0.0},
        {"dynamic_viscosity": -1e-3},
        {"yield_stress": -1.0},
        {"h_viscous_min": 0.0},
    ],
)
def test_solver_config_rejects_invalid_viscous_parameters(kwargs: dict[str, float | str]) -> None:
    with pytest.raises(ValueError):
        SolverConfig(nx=2, ny=2, dx=1.0, dy=1.0, **kwargs)
