"""Diagnostic tools for factorial designs: Main effects, interactions, and half-normal probability plots."""

from __future__ import annotations

import itertools
from typing import Sequence, Union

import numpy as np
import pandas as pd
from scipy import stats

from doe_toolkit.base import DesignMatrix


def calculate_effects(design: DesignMatrix, y: Sequence[float] | np.ndarray) -> pd.DataFrame:
    """Calculate main effects, 2-factor interactions, and sum of squares for 2-level designs."""
    X = design.coded_matrix
    y_arr = np.asarray(y, dtype=np.float64).ravel()
    n_runs = len(y_arr)

    if len(X) != n_runs:
        raise ValueError(f"Design rows ({len(X)}) must match response length ({n_runs}).")

    factors = design.factor_names
    k = len(factors)
    ss_total = float(np.sum((y_arr - np.mean(y_arr)) ** 2))

    rows = []

    # 1. Main effects
    for idx, factor in enumerate(factors):
        x_col = X[:, idx]
        # Contrast sum
        contrast = float(np.sum(x_col * y_arr))
        # Number of plus/minus runs is n_runs / 2 for orthogonal 2-level designs
        n_half = n_runs / 2.0
        effect = contrast / n_half
        coef = effect / 2.0
        ss = (contrast ** 2) / n_runs
        pct_contrib = (ss / ss_total * 100.0) if ss_total > 0 else 0.0

        rows.append(
            {
                "Term": factor,
                "Type": "Main Effect",
                "Effect": effect,
                "Coef": coef,
                "SS": ss,
                "Percent_Contribution": pct_contrib,
            }
        )

    # 2. Two-factor interaction effects
    for idx1, idx2 in itertools.combinations(range(k), 2):
        f1, f2 = factors[idx1], factors[idx2]
        x_inter = X[:, idx1] * X[:, idx2]
        contrast = float(np.sum(x_inter * y_arr))
        n_half = n_runs / 2.0
        effect = contrast / n_half
        coef = effect / 2.0
        ss = (contrast ** 2) / n_runs
        pct_contrib = (ss / ss_total * 100.0) if ss_total > 0 else 0.0

        rows.append(
            {
                "Term": f"{f1}:{f2}",
                "Type": "2-Factor Interaction",
                "Effect": effect,
                "Coef": coef,
                "SS": ss,
                "Percent_Contribution": pct_contrib,
            }
        )

    df_effects = pd.DataFrame(rows).sort_values("SS", ascending=False).reset_index(drop=True)
    return df_effects


def half_normal_plot_coordinates(effects: Union[pd.DataFrame, dict[str, float]]) -> pd.DataFrame:
    """Compute coordinates for Daniel's Half-Normal Probability Plot to identify active effects without replicates."""
    if isinstance(effects, pd.DataFrame):
        terms = effects["Term"].tolist()
        vals = effects["Effect"].tolist()
    else:
        terms = list(effects.keys())
        vals = list(effects.values())

    abs_vals = [abs(v) for v in vals]
    m = len(terms)

    # Sort in ascending order of absolute effect
    sorted_indices = np.argsort(abs_vals)
    sorted_terms = [terms[i] for i in sorted_indices]
    sorted_abs = [abs_vals[i] for i in sorted_indices]
    sorted_raw = [vals[i] for i in sorted_indices]

    # Half-normal probabilities: p_i = (i - 0.5) / m
    i_ranks = np.arange(1, m + 1)
    probs = (i_ranks - 0.5) / m
    # Half-normal quantiles: z_i = Phi^-1(0.5 + 0.5 * p_i)
    z_quantiles = stats.norm.ppf(0.5 + 0.5 * probs)

    return pd.DataFrame(
        {
            "Rank": i_ranks,
            "Term": sorted_terms,
            "Raw_Effect": sorted_raw,
            "Absolute_Effect": sorted_abs,
            "Probability": probs,
            "Half_Normal_Quantile": z_quantiles,
        }
    )
