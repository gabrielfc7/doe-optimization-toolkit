"""Base classes and container types for experimental designs."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd


@dataclass
class DesignMatrix:
    """Container for experimental design matrices in coded and natural units."""

    coded_matrix: np.ndarray
    factor_names: list[str]
    factor_levels: Optional[dict[str, tuple[float, float]]] = None
    design_type: str = "custom"
    extra_info: dict[str, Any] = field(default_factory=dict)

    @property
    def n_runs(self) -> int:
        return len(self.coded_matrix)

    @property
    def n_factors(self) -> int:
        return len(self.factor_names)

    def to_dataframe(self) -> pd.DataFrame:
        """Return design matrix in coded units (-1, 0, 1)."""
        return pd.DataFrame(self.coded_matrix, columns=self.factor_names)

    def to_natural_dataframe(self) -> pd.DataFrame:
        """Convert coded matrix to engineering / natural units based on factor_levels."""
        if self.factor_levels is None:
            return self.to_dataframe()

        df_nat = pd.DataFrame(index=range(self.n_runs))
        for factor in self.factor_names:
            low, high = self.factor_levels.get(factor, (-1.0, 1.0))
            center = (high + low) / 2.0
            half_range = (high - low) / 2.0
            coded_col = self.coded_matrix[:, self.factor_names.index(factor)]
            df_nat[factor] = center + coded_col * half_range

        return df_nat

    def __str__(self) -> str:
        return self.to_dataframe().to_string()
