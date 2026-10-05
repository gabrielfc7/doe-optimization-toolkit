# Mathematical Foundations: Effects & Half-Normal Plots

This document details effect calculations and sparsity-of-effects screening implemented in `doe_toolkit.diagnostics`.

---

## 1. Contrast & Effect Estimation

For an orthogonal two-level design with $N$ runs and response vector $y \in \mathbb{R}^N$:

### 1.1 Main Effects

The main effect of factor $i$ is the average change in response produced by moving factor $i$ from its low level ($-1$) to its high level ($+1$):

$$ME_i = \bar{y}_{i+} - \bar{y}_{i-} = \frac{1}{N/2} \sum_{j=1}^N x_{ji} y_j = \frac{\text{Contrast}_i}{N/2}$$

Regression coefficient $\beta_i = \frac{ME_i}{2}$.

### 1.2 Two-Factor Interactions

$$INT_{im} = \frac{1}{N/2} \sum_{j=1}^N (x_{ji} x_{jm}) y_j = \frac{\text{Contrast}_{im}}{N/2}$$

### 1.3 Sum of Squares & Variance Contribution

$$SS_E = \frac{(\text{Contrast}_E)^2}{N} = \frac{N}{4} (\text{Effect}_E)^2$$

Given total corrected sum of squares $SS_{\text{total}} = \sum (y_j - \bar{y})^2$:

$$\text{Contribution}_E = \frac{SS_E}{SS_{\text{total}}} \times 100\%$$

---

## 2. Daniel's Half-Normal Probability Plot

In unreplicated factorial experiments ($N$ runs with $N - 1$ estimated effects), there are zero degrees of freedom left to estimate experimental error variance $\sigma^2$.

By the **sparsity-of-effects principle**, most effects in a screening experiment are inactive noise:

$$\text{Effect}_i \sim \mathcal{N}(0, \sigma_E^2)$$

### 2.1 Theoretical Half-Normal Distribution

If $Z \sim \mathcal{N}(0, \sigma^2)$, then $|Z|$ follows a half-normal distribution with CDF:

$$F_{|Z|}(z) = 2 \Phi\left(\frac{z}{\sigma}\right) - 1, \quad z \ge 0$$

### 2.2 Quantile Mapping

1. Rank the $m$ absolute effects in ascending order: $|E|_{(1)} \le |E|_{(2)} \le \dots \le |E|_{(m)}$.
2. Empirical plotting positions:
   $$p_i = \frac{i - 0.5}{m}$$
3. Standard half-normal quantiles:
   $$z_i = \Phi^{-1}\left( 0.5 + 0.5 p_i \right)$$

Inactive effects fall along a straight line passing through the origin. Active (statistically significant) effects deviate substantially above the reference line.
