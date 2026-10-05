"""Tests for factorial and fractional factorial design generation."""

import numpy as np
import pytest

from doe_toolkit.factorial import (
    fractional_factorial_2k_p,
    full_factorial_2k,
    full_factorial_3k,
)


def test_full_factorial_2k_orthogonality():
    factors = ["A", "B", "C"]
    design = full_factorial_2k(factors)
    X = design.coded_matrix

    assert X.shape == (8, 3)
    # Check orthogonality: X^T * X = N * I
    xt_x = X.T @ X
    expected = 8.0 * np.eye(3)
    np.testing.assert_allclose(xt_x, expected, atol=1e-10)


def test_full_factorial_3k_levels():
    factors = ["A", "B"]
    design = full_factorial_3k(factors)
    X = design.coded_matrix

    assert X.shape == (9, 2)
    assert set(np.unique(X)) == {-1.0, 0.0, 1.0}


def test_fractional_factorial_resolution_and_aliases():
    # 2^(4-1) design with generator D = A*B*C -> Defining relation I = ABCD (Resolution IV)
    factors = ["A", "B", "C", "D"]
    generators = {"D": "A*B*C"}
    design, alias_map, resolution = fractional_factorial_2k_p(factors, generators)

    assert design.n_runs == 8
    assert "Resolution IV" in resolution
    assert design.extra_info["defining_relation"] == ["ABCD"]

    # Verify alias relationships:
    # A is aliased with BCD
    assert "BCD" in alias_map["A"]
    # AB is aliased with CD
    assert "CD" in alias_map["AB"]
    # BC is aliased with AD
    assert "AD" in alias_map["BC"]


def test_fractional_factorial_resolution_iii():
    # 2^(3-1) design with C = A*B -> Defining relation I = ABC (Resolution III)
    factors = ["A", "B", "C"]
    generators = {"C": "A*B"}
    design, alias_map, resolution = fractional_factorial_2k_p(factors, generators)

    assert design.n_runs == 4
    assert "Resolution III" in resolution
    assert "BC" in alias_map["A"]
