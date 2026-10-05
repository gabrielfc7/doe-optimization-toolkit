# Mathematical Foundations: Response Surface Methodology (RSM)

This document formalizes second-order experimental design and stationary point optimization implemented in `doe_toolkit.rsm`.

---

## 1. Central Composite Design (CCD)

A Central Composite Design augments a $2^k$ factorial cube with $2k$ axial (star) points and $n_c$ center points.

### 1.1 Axial Distance ($\alpha$) Selection

* **Rotatability:** A design is rotatable if the prediction variance $\text{Var}(\hat{y}(x))$ depends only on the distance from the design center $\rho = \|x\|$, and not on the direction. Rotatability requires:
  $$\alpha = (N_F)^{1/4} = (2^k)^{1/4}$$
* **Spherical Design:** All factorial and axial points lie on a sphere of radius $\sqrt{k}$:
  $$\alpha = \sqrt{k}$$
* **Face-Centered Design:** Axial points lie on the center of the faces of the factorial space ($\alpha = 1.0$), requiring only 3 levels per factor.

---

## 2. Second-Order Polynomial Model

$$\hat{y}(x) = \beta_0 + \sum_{i=1}^k \beta_i x_i + \sum_{i=1}^k \beta_{ii} x_i^2 + \sum_{i < j} \beta_{ij} x_i x_j$$

In vector-matrix notation:

$$\hat{y}(x) = \beta_0 + x^T b + x^T B x$$

where:
* $b = [\beta_1, \beta_2, \dots, \beta_k]^T \in \mathbb{R}^k$ (linear gradient vector).
* $B \in \mathbb{R}^{k \times k}$ is the symmetric second-order Hessian matrix:
  $$B_{ii} = \beta_{ii}, \quad B_{ij} = B_{ji} = \frac{1}{2} \beta_{ij} \quad (i \ne j)$$

---

## 3. Stationary Point Optimization & Canonical Analysis

Taking the derivative with respect to $x$ and setting to zero:

$$\frac{\partial \hat{y}}{\partial x} = b + 2 B x = \mathbf{0}$$

### 3.1 Stationary Point Location

$$x_s = -\frac{1}{2} B^{-1} b$$

Predicted response at the stationary point:

$$\hat{y}_s = \beta_0 + \frac{1}{2} x_s^T b$$

### 3.2 Canonical Classification via Hessian Eigenvalues

Let $B = M \Lambda M^T$ be the spectral decomposition of $B$ with eigenvalues $\lambda_1, \dots, \lambda_k$:

* **Local Maximum:** All $\lambda_i < 0$ (negative definite $B$). The stationary point is a peak.
* **Local Minimum:** All $\lambda_i > 0$ (positive definite $B$). The stationary point is a valley.
* **Saddle Point:** Eigenvalues have mixed signs (both positive and negative). The response surface is a saddle (minimax surface).
