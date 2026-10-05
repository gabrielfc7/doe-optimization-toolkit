"""Factorial design generators: 2^k, 3^k, and 2^(k-p) fractional factorials with alias structures."""

from __future__ import annotations

import itertools
from typing import Optional

import numpy as np

from doe_toolkit.base import DesignMatrix


def full_factorial_2k(
    factors: list[str],
    center_points: int = 0,
    factor_levels: Optional[dict[str, tuple[float, float]]] = None,
) -> DesignMatrix:
    """Generate a standard 2^k two-level full factorial design matrix in coded units (-1, +1)."""
    k = len(factors)
    if k < 1:
        raise ValueError("Must provide at least 1 factor.")

    # Generate all 2^k combinations in standard order: [-1, +1]
    levels = list(itertools.product([-1.0, 1.0], repeat=k))
    matrix = np.array(levels, dtype=np.float64)

    if center_points > 0:
        centers = np.zeros((center_points, k), dtype=np.float64)
        matrix = np.vstack([matrix, centers])

    return DesignMatrix(
        coded_matrix=matrix,
        factor_names=factors,
        factor_levels=factor_levels,
        design_type=f"2^{k} Full Factorial",
        extra_info={"center_points": center_points},
    )


def full_factorial_3k(
    factors: list[str],
    factor_levels: Optional[dict[str, tuple[float, float]]] = None,
) -> DesignMatrix:
    """Generate a 3^k three-level full factorial design matrix in coded units (-1, 0, +1)."""
    k = len(factors)
    levels = list(itertools.product([-1.0, 0.0, 1.0], repeat=k))
    matrix = np.array(levels, dtype=np.float64)

    return DesignMatrix(
        coded_matrix=matrix,
        factor_names=factors,
        factor_levels=factor_levels,
        design_type=f"3^{k} Full Factorial",
    )


def _multiply_words(w1: str, w2: str) -> str:
    """Multiply two factor words using group identity A * A = I."""
    chars1 = set(w1.replace("*", ""))
    chars2 = set(w2.replace("*", ""))
    sym_diff = sorted(chars1 ^ chars2)
    return "".join(sym_diff) if sym_diff else "I"


def fractional_factorial_2k_p(
    factors: list[str],
    generators: dict[str, str],
    factor_levels: Optional[dict[str, tuple[float, float]]] = None,
) -> tuple[DesignMatrix, dict[str, list[str]], str]:
    """Generate a 2^(k-p) fractional factorial design with alias structure and resolution.

    Parameters
    ----------
    factors : list of str
        List of all factors (e.g. ['A', 'B', 'C', 'D']).
    generators : dict of str -> str
        Generators for added factors (e.g. {'D': 'A*B*C'} or {'D': 'A*B', 'E': 'A*C'}).
    """
    p = len(generators)
    k = len(factors)
    base_factors = [f for f in factors if f not in generators]
    k_base = len(base_factors)

    if k_base + p != k:
        raise ValueError("Total factors must equal base factors plus generator factors.")

    # 1. Generate full factorial on base factors
    base_design = full_factorial_2k(base_factors).coded_matrix
    cols = {f: base_design[:, idx] for idx, f in enumerate(base_factors)}

    # 2. Compute generated factors as products
    for gen_factor, gen_formula in generators.items():
        terms = [t.strip() for t in gen_formula.split("*")]
        prod_col = np.ones(len(base_design), dtype=np.float64)
        for t in terms:
            if t not in cols:
                raise ValueError(f"Term '{t}' in generator '{gen_formula}' is not in available factors.")
            prod_col *= cols[t]
        cols[gen_factor] = prod_col

    # Construct final matrix in order of factors
    full_matrix = np.column_stack([cols[f] for f in factors])

    # 3. Defining relation: for each generator D = ABC -> Word = ABCD
    base_words = []
    for gen_factor, gen_formula in generators.items():
        word = _multiply_words(gen_factor, gen_formula)
        base_words.append(word)

    # All combinations of defining words to find complete defining relation
    all_defining_words = set(base_words)
    for r in range(2, len(base_words) + 1):
        for combo in itertools.combinations(base_words, r):
            curr = combo[0]
            for w in combo[1:]:
                curr = _multiply_words(curr, w)
            if curr != "I":
                all_defining_words.add(curr)

    sorted_words = sorted(all_defining_words, key=lambda w: (len(w), w))
    resolution_num = min(len(w) for w in sorted_words) if sorted_words else k
    roman_numerals = {3: "III", 4: "IV", 5: "V", 6: "VI", 7: "VII"}
    resolution_str = f"Resolution {roman_numerals.get(resolution_num, str(resolution_num))}"

    # 4. Alias structure for main effects and 2-factor interactions
    alias_map: dict[str, list[str]] = {}
    candidate_effects = list(factors) + ["".join(sorted(pair)) for pair in itertools.combinations(factors, 2)]

    for effect in candidate_effects:
        aliases = []
        for def_word in sorted_words:
            alias = _multiply_words(effect, def_word)
            if alias != "I" and alias != effect:
                aliases.append(alias)
        alias_map[effect] = sorted(aliases, key=lambda w: (len(w), w))

    design = DesignMatrix(
        coded_matrix=full_matrix,
        factor_names=factors,
        factor_levels=factor_levels,
        design_type=f"2^({k}-{p}) Fractional Factorial ({resolution_str})",
        extra_info={"defining_relation": sorted_words, "resolution": resolution_str},
    )

    return design, alias_map, resolution_str
