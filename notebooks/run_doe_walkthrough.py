"""Demonstration and walkthrough script for doe-optimization-toolkit."""

import sys
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np

from doe_toolkit.factorial import full_factorial_2k, fractional_factorial_2k_p
from doe_toolkit.diagnostics import calculate_effects, half_normal_plot_coordinates
from doe_toolkit.rsm import central_composite_design, QuadraticResponseSurface


def main():
    print("=" * 80)
    print("           DOE-OPTIMIZATION-TOOLKIT DEMONSTRATION & VERIFICATION")
    print("=" * 80)

    # 1. Full Factorial 2^k Design & Engineering Units
    factors = ["Temperature", "Pressure", "Catalyst"]
    levels = {"Temperature": (120.0, 180.0), "Pressure": (2.0, 6.0), "Catalyst": (0.5, 2.5)}
    full_design = full_factorial_2k(factors, center_points=2, factor_levels=levels)

    print("\n### 1. Full 2^3 Factorial Design Matrix in Natural Engineering Units:")
    print(full_design.to_natural_dataframe())

    # 2. Fractional Factorial 2^(4-1) with Defining Relation & Alias Structure
    factors_4 = ["A", "B", "C", "D"]
    generators = {"D": "A*B*C"}
    frac_design, alias_map, resolution = fractional_factorial_2k_p(factors_4, generators)

    print(f"\n### 2. Fractional Factorial Design ({resolution}):")
    print(f"  Defining Relation: I = {' = '.join(frac_design.extra_info['defining_relation'])}")
    print("  Confounded Alias Structure:")
    for effect, aliases in list(alias_map.items())[:5]:
        print(f"    {effect} = {' + '.join(aliases)}")

    # 3. Factorial Diagnostic Effects & Half-Normal Plot Coordinates
    X = full_design.coded_matrix[:8]
    # Simulated response: Yield = 60 + 10*Temp - 4*Press + 14*Cat + 5*(Temp*Press)
    y = 60.0 + 10.0 * X[:, 0] - 4.0 * X[:, 1] + 14.0 * X[:, 2] + 5.0 * (X[:, 0] * X[:, 1])

    sub_design = full_factorial_2k(factors)
    effects_df = calculate_effects(sub_design, y)
    print("\n### 3. Factorial Main & Interaction Effects Summary:")
    print(effects_df.to_string(index=False))

    hn_df = half_normal_plot_coordinates(effects_df)
    print("\n### 4. Daniel's Half-Normal Plot Coordinates (Top Active Factors):")
    print(hn_df.tail(3).to_string(index=False))

    # 4. Central Composite Design & Response Surface Optimization
    rsm_factors = ["Temp", "Conc"]
    rsm_levels = {"Temp": (140.0, 180.0), "Conc": (10.0, 30.0)}
    ccd = central_composite_design(rsm_factors, alpha="rotatable", center_points=5, factor_levels=rsm_levels)

    # True surface with maximum at coded point (0.3, -0.4)
    # y = 92.0 - 6*(x1 - 0.3)^2 - 8*(x2 + 0.4)^2
    X_ccd = ccd.coded_matrix
    y_rsm = 92.0 - 6.0 * (X_ccd[:, 0] - 0.3) ** 2 - 8.0 * (X_ccd[:, 1] + 0.4) ** 2

    surface_model = QuadraticResponseSurface()
    surface_model.fit(ccd, y_rsm)
    opt_result = surface_model.optimize()

    print("\n### 5. Response Surface Optimization & Canonical Analysis:")
    print(f"  Surface R^2: {surface_model.r2_:.4f}")
    print(f"  Stationary Point (Coded Units):   Temp = {opt_result.point_coded[0]:.4f}, Conc = {opt_result.point_coded[1]:.4f}")
    if opt_result.point_natural:
        print(f"  Stationary Point (Natural Units): Temp = {opt_result.point_natural['Temp']:.2f} °C, Conc = {opt_result.point_natural['Conc']:.2f}%")
    print(f"  Predicted Optimal Response:       {opt_result.predicted_response:.2f}%")
    print(f"  Nature of Stationary Point:       {opt_result.nature} (Eigenvalues: {np.round(opt_result.eigenvalues, 3)})")
    print("=" * 80)


if __name__ == "__main__":
    main()
