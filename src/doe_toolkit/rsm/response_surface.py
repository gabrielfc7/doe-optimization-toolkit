"""Response Surface Methodology (RSM): CCD, Box-Behnken, and Quadratic Optimization."""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Literal, Optional, Sequence, Union

import numpy as np
import pandas as pd

from doe_toolkit.base import DesignMatrix


def central_composite_design(
    factors: list[str],
    alpha: Union[Literal["rotatable", "spherical", "face_centered"], float] = "rotatable",
    center_points: int = 4,
    factor_levels: Optional[dict[str, tuple[float, float]]] = None,
) -> DesignMatrix:
    """Generate a Central Composite Design (CCD) for response surface modeling.

    Parameters
    ----------
    factors : list of str
        Names of factors.
    alpha : {'rotatable', 'spherical', 'face_centered'} or float, default='rotatable'
        Distance of the axial (star) points from the design center.
    center_points : int, default=4
        Number of center points at (0, 0, ..., 0).
    """
    k = len(factors)
    if k < 2:
        raise ValueError("CCD requires at least 2 factors.")

    # 1. Factorial points: 2^k
    factorial_points = np.array(list(itertools.product([-1.0, 1.0], repeat=k)), dtype=np.float64)
    n_fact = len(factorial_points)

    # 2. Compute alpha value
    if alpha == "rotatable":
        alpha_val = float((n_fact) ** 0.25)
    elif alpha == "spherical":
        alpha_val = float(np.sqrt(k))
    elif alpha == "face_centered":
        alpha_val = 1.0
    elif isinstance(alpha, (int, float)):
        alpha_val = float(alpha)
    else:
        raise ValueError(f"Unknown alpha parameter '{alpha}'.")

    # 3. Axial / star points: 2k
    axial_points = np.zeros((2 * k, k), dtype=np.float64)
    for i in range(k):
        axial_points[2 * i, i] = -alpha_val
        axial_points[2 * i + 1, i] = alpha_val

    # 4. Center points
    center_pts = np.zeros((center_points, k), dtype=np.float64)

    full_matrix = np.vstack([factorial_points, axial_points, center_pts])

    return DesignMatrix(
        coded_matrix=full_matrix,
        factor_names=factors,
        factor_levels=factor_levels,
        design_type=f"Central Composite Design (alpha={alpha_val:.3f})",
        extra_info={"alpha": alpha_val, "n_factorial": n_fact, "n_axial": 2 * k, "n_center": center_points},
    )


def box_behnken_design(
    factors: list[str],
    center_points: int = 3,
    factor_levels: Optional[dict[str, tuple[float, float]]] = None,
) -> DesignMatrix:
    """Generate a Box-Behnken Design (BBD) with 3 levels per factor (-1, 0, +1)."""
    k = len(factors)
    if k < 3:
        raise ValueError("Box-Behnken design requires at least 3 factors.")

    blocks = []
    # For every pair of factors, generate a 2^2 full factorial while holding other factors at 0
    for i, j in itertools.combinations(range(k), 2):
        for s1, s2 in itertools.product([-1.0, 1.0], repeat=2):
            row = np.zeros(k, dtype=np.float64)
            row[i] = s1
            row[j] = s2
            blocks.append(row)

    edge_matrix = np.array(blocks, dtype=np.float64)
    center_pts = np.zeros((center_points, k), dtype=np.float64)
    full_matrix = np.vstack([edge_matrix, center_pts])

    return DesignMatrix(
        coded_matrix=full_matrix,
        factor_names=factors,
        factor_levels=factor_levels,
        design_type="Box-Behnken Design",
        extra_info={"n_edge": len(edge_matrix), "n_center": center_points},
    )


@dataclass
class StationaryPointResult:
    """Container for response surface stationary point and canonical analysis."""

    point_coded: np.ndarray
    point_natural: Optional[dict[str, float]]
    predicted_response: float
    nature: str  # 'Maximum', 'Minimum', 'Saddle Point'
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray


class QuadraticResponseSurface:
    """Quadratic response surface model with stationary point optimization and canonical analysis."""

    def __init__(self):
        self.factor_names: list[str] = []
        self.factor_levels: Optional[dict[str, tuple[float, float]]] = None
        self.coef_: Optional[np.ndarray] = None
        self.intercept_: float = 0.0
        self.r2_: float = 0.0
        self.feature_names_: list[str] = []
        self.b_linear_: Optional[np.ndarray] = None
        self.B_hessian_: Optional[np.ndarray] = None

    def _build_features(self, X: np.ndarray) -> np.ndarray:
        n_samples, k = X.shape
        cols = [np.ones(n_samples, dtype=np.float64)]

        # Linear terms
        for i in range(k):
            cols.append(X[:, i])

        # Quadratic terms
        for i in range(k):
            cols.append(X[:, i] ** 2)

        # Cross interaction terms
        for i, j in itertools.combinations(range(k), 2):
            cols.append(X[:, i] * X[:, j])

        return np.column_stack(cols)

    def fit(self, design: DesignMatrix, y: Sequence[float] | np.ndarray) -> QuadraticResponseSurface:
        """Fit a second-order polynomial response surface to experimental runs."""
        X_coded = design.coded_matrix
        self.factor_names = design.factor_names
        self.factor_levels = design.factor_levels
        y_arr = np.asarray(y, dtype=np.float64).ravel()

        k = len(self.factor_names)
        self.feature_names_ = ["Intercept"] + self.factor_names + [f"{f}^2" for f in self.factor_names]
        self.feature_names_ += [f"{self.factor_names[i]}:{self.factor_names[j]}" for i, j in itertools.combinations(range(k), 2)]

        X_poly = self._build_features(X_coded)

        # Ordinary Least Squares: beta = (X^T X)^-1 X^T y
        beta, _, _, _ = np.linalg.lstsq(X_poly, y_arr, rcond=None)

        self.coef_ = beta
        self.intercept_ = float(beta[0])

        # Linear gradient vector b
        self.b_linear_ = beta[1 : 1 + k]

        # Symmetric second-order Hessian matrix B
        self.B_hessian_ = np.zeros((k, k), dtype=np.float64)
        for i in range(k):
            self.B_hessian_[i, i] = beta[1 + k + i]  # quadratic coefficient beta_ii

        inter_idx = 1 + 2 * k
        for i, j in itertools.combinations(range(k), 2):
            val = beta[inter_idx] / 2.0  # off-diagonal B_ij = 0.5 * beta_ij
            self.B_hessian_[i, j] = val
            self.B_hessian_[j, i] = val
            inter_idx += 1

        # R^2 calculation
        y_pred = X_poly @ beta
        ss_res = np.sum((y_arr - y_pred) ** 2)
        ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)
        self.r2_ = float(1.0 - (ss_res / ss_tot)) if ss_tot > 0 else 0.0

        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict response for points in coded units."""
        X_poly = self._build_features(X)
        return X_poly @ self.coef_

    def optimize(self) -> StationaryPointResult:
        """Find the stationary point x_s = -0.5 * B^-1 * b and classify its nature."""
        if self.B_hessian_ is None or self.b_linear_ is None:
            raise RuntimeError("Model must be fitted before optimization.")

        # Solve B * x_s = -0.5 * b
        try:
            x_s = np.linalg.solve(self.B_hessian_, -0.5 * self.b_linear_)
        except np.linalg.LinAlgError:
            x_s = np.linalg.pinv(self.B_hessian_) @ (-0.5 * self.b_linear_)

        # Predicted response at stationary point: y_s = beta_0 + 0.5 * x_s^T * b
        y_s = float(self.intercept_ + 0.5 * np.dot(x_s, self.b_linear_))

        # Canonical analysis via eigenvalues of B
        eigenvalues, eigenvectors = np.linalg.eigh(self.B_hessian_)

        tol = 1e-6
        if np.all(eigenvalues < -tol):
            nature = "Maximum"
        elif np.all(eigenvalues > tol):
            nature = "Minimum"
        else:
            nature = "Saddle Point"

        # Convert to natural units if available
        nat_point = None
        if self.factor_levels:
            nat_point = {}
            for idx, factor in enumerate(self.factor_names):
                low, high = self.factor_levels.get(factor, (-1.0, 1.0))
                center = (high + low) / 2.0
                half_range = (high - low) / 2.0
                nat_point[factor] = float(center + x_s[idx] * half_range)

        return StationaryPointResult(
            point_coded=x_s,
            point_natural=nat_point,
            predicted_response=y_s,
            nature=nature,
            eigenvalues=eigenvalues,
            eigenvectors=eigenvectors,
        )
