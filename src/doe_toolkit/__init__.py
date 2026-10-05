"""doe-optimization-toolkit: Design of Experiments & Response Surface Methodology."""

__version__ = "0.1.0"

from doe_toolkit.base import DesignMatrix
from doe_toolkit.diagnostics import (
    calculate_effects,
    half_normal_plot_coordinates,
)
from doe_toolkit.factorial import (
    fractional_factorial_2k_p,
    full_factorial_2k,
    full_factorial_3k,
)
from doe_toolkit.rsm import (
    QuadraticResponseSurface,
    StationaryPointResult,
    box_behnken_design,
    central_composite_design,
)

__all__ = [
    "DesignMatrix",
    "full_factorial_2k",
    "full_factorial_3k",
    "fractional_factorial_2k_p",
    "calculate_effects",
    "half_normal_plot_coordinates",
    "central_composite_design",
    "box_behnken_design",
    "QuadraticResponseSurface",
    "StationaryPointResult",
]
