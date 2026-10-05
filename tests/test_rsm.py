"""Tests for Response Surface Methodology (CCD, Box-Behnken, and Quadratic Surface)."""

import numpy as np
import pytest

from doe_toolkit.rsm import (
    QuadraticResponseSurface,
    box_behnken_design,
    central_composite_design,
)


def test_central_composite_design_structure():
    factors = ["X1", "X2", "X3"]
    # 2^3 = 8 factorial + 2*3 = 6 axial + 4 center = 18 runs
    design = central_composite_design(factors, alpha="rotatable", center_points=4)

    assert design.n_runs == 18
    assert design.n_factors == 3
    # Rotatable alpha for k=3 is (8)^0.25 = 1.6818
    assert pytest.approx(design.extra_info["alpha"], abs=1e-3) == 1.682


def test_box_behnken_design_structure():
    factors = ["X1", "X2", "X3"]
    # For k=3, 12 edge runs + 3 center runs = 15 runs
    design = box_behnken_design(factors, center_points=3)

    assert design.n_runs == 15
    assert set(np.unique(design.coded_matrix)) == {-1.0, 0.0, 1.0}


def test_quadratic_response_surface_optimization(rsm_dome_surface):
    design, y = rsm_dome_surface
    model = QuadraticResponseSurface()
    model.fit(design, y)

    # Perfect quadratic surface -> R^2 should be ~ 1.0
    assert model.r2_ > 0.999

    # Optimize and classify stationary point
    res = model.optimize()
    assert res.nature == "Maximum"
    # Ground truth peak was (0.2, -0.3)
    np.testing.assert_allclose(res.point_coded, [0.2, -0.3], atol=1e-4)
    # Ground truth peak response was 80.0
    assert pytest.approx(res.predicted_response, abs=1e-3) == 80.0
    # Eigenvalues of Hessian should both be negative for a local maximum
    assert np.all(res.eigenvalues < 0)
