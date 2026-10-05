"""Test fixtures for doe-optimization-toolkit."""

import numpy as np
import pytest

from doe_toolkit.factorial import full_factorial_2k
from doe_toolkit.rsm import central_composite_design


@pytest.fixture
def chemical_yield_3k():
    """3-factor experiment with temperature, pressure, and catalyst."""
    factors = ["Temp", "Press", "Cat"]
    design = full_factorial_2k(factors)
    X = design.coded_matrix

    # Synthetic response: true model y = 50 + 8*Temp - 4*Press + 12*Cat + 6*(Temp*Press) + noise
    rng = np.random.RandomState(42)
    y = 50.0 + 8.0 * X[:, 0] - 4.0 * X[:, 1] + 12.0 * X[:, 2] + 6.0 * (X[:, 0] * X[:, 1]) + rng.normal(0, 0.5, size=len(X))
    return design, y


@pytest.fixture
def rsm_dome_surface():
    """2-factor quadratic response surface with known local maximum at (0.2, -0.3)."""
    factors = ["X1", "X2"]
    design = central_composite_design(factors, alpha="rotatable", center_points=4)
    X = design.coded_matrix

    # Quadratic response: y = 80 - 5*(x1 - 0.2)^2 - 8*(x2 + 0.3)^2
    # Expanding: 80 - 5*(x1^2 - 0.4*x1 + 0.04) - 8*(x2^2 + 0.6*x2 + 0.09)
    # y = 79.08 + 2.0*x1 - 4.8*x2 - 5*x1^2 - 8*x2^2
    y = 79.08 + 2.0 * X[:, 0] - 4.8 * X[:, 1] - 5.0 * (X[:, 0] ** 2) - 8.0 * (X[:, 1] ** 2)
    return design, y
