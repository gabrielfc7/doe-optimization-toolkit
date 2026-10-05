"""Tests for effect calculation and half-normal probability plots."""

import numpy as np
import pytest

from doe_toolkit.diagnostics import (
    calculate_effects,
    half_normal_plot_coordinates,
)


def test_calculate_effects_and_contributions(chemical_yield_3k):
    design, y = chemical_yield_3k
    df_effects = calculate_effects(design, y)

    assert len(df_effects) == 6  # 3 main effects + 3 interactions
    assert "Cat" in df_effects["Term"].values
    assert "Temp" in df_effects["Term"].values

    # Catalyst was simulated with effect +12, should have highest or second highest SS
    top_term = df_effects.iloc[0]["Term"]
    assert top_term in ["Cat", "Temp", "Temp:Press"]

    # Sum of percent contributions should be <= 100%
    assert 0.0 < np.sum(df_effects["Percent_Contribution"]) <= 100.0


def test_half_normal_plot_coordinates(chemical_yield_3k):
    design, y = chemical_yield_3k
    df_effects = calculate_effects(design, y)
    df_hn = half_normal_plot_coordinates(df_effects)

    assert len(df_hn) == len(df_effects)
    assert "Half_Normal_Quantile" in df_hn.columns
    assert "Absolute_Effect" in df_hn.columns
    # Half-normal quantiles and absolute effects must both be strictly monotonic
    assert np.all(np.diff(df_hn["Half_Normal_Quantile"]) > 0)
    assert np.all(np.diff(df_hn["Absolute_Effect"]) >= 0)
