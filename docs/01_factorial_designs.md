# Mathematical Foundations: Full and Fractional Factorial Designs

This document formalizes design generation, confounding mechanics, and alias structures implemented in `doe_toolkit.factorial`.

---

## 1. Two-Level Full Factorial Designs ($2^k$)

A $2^k$ factorial design evaluates all $2^k$ combinations of $k$ factors at two levels (coded as $-1$ and $+1$).

### 1.1 Orthogonality Property

Let $X \in \{-1, +1\}^{N \times k}$ denote the coded design matrix with $N = 2^k$. Every column $x_i$ satisfies:

$$\sum_{j=1}^N x_{ji} = 0, \quad \sum_{j=1}^N x_{ji}^2 = N$$

For any distinct pair of factors $i \ne m$:

$$\sum_{j=1}^N x_{ji} x_{jm} = 0 \iff X^T X = N \cdot I_k$$

Orthogonality guarantees that every main effect and interaction estimate is statistically uncorrelated with zero covariance: $\text{Cov}(\hat{\beta}_i, \hat{\beta}_m) = 0$.

---

## 2. Fractional Factorial Designs ($2^{k-p}$)

When running all $2^k$ combinations is cost-prohibitive, a $2^{k-p}$ fraction evaluates $\frac{1}{2^p}$ of the full design.

### 2.1 Generators & Defining Relation

We select $k - p$ base factors and assign the remaining $p$ factors to interaction columns of the base factors.

For example, in a $2^{4-1}$ design ($k=4, p=1$):
* Base factors: $A, B, C$ ($2^3 = 8$ runs)
* Generator: $D = ABC$

Multiplying both sides by $D$ using modulo-2 exponentiation ($D^2 \equiv I$):

$$I = ABCD$$

The word $ABCD$ is the **defining relation**.

### 2.2 Design Resolution

The resolution of a fractional factorial design corresponds to the length of the shortest word in the defining relation:

* **Resolution III:** No main effects are aliased with other main effects, but main effects are aliased with 2-factor interactions ($A = BC$).
* **Resolution IV:** No main effect is aliased with any 2-factor interaction, but 2-factor interactions are aliased with each other ($AB = CD$).
* **Resolution V:** No main effect or 2-factor interaction is aliased with any other main effect or 2-factor interaction. Two-factor interactions are aliased only with 3-factor interactions.

### 2.3 Alias Structure

Multiplying any candidate effect $E$ by each word in the defining relation identifies all confounded effects:

$$E \cdot I = E \cdot (ABCD) \implies E = E \cdot ABCD$$

For main effect $A$:

$$A \cdot (ABCD) = A^2 BCD = I \cdot BCD = BCD \implies A \equiv BCD$$

For two-factor interaction $AB$:

$$AB \cdot (ABCD) = A^2 B^2 CD = CD \implies AB \equiv CD$$
