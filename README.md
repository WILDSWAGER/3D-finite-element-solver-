# 3D Linear Static Finite Element Analysis (FEA) for Nodal Displacement

This repository contains a standalone, parallelized MATLAB script for performing 3D linear static structural finite element analysis (FEA). The algorithm imports a customized CAD geometry from an STL file, discretizes the geometry into linear tetrahedral elements, assembles a global stiffness matrix, and evaluates the corresponding 3D displacement vector fields along with boundary reaction forces under localized compressive loading.

## 🛠 Features

* **STL Geometry Ingestion:** Utilizes MATLAB's Partial Differential Equation (PDE) Toolbox to import and visualize 3D `.stl` continuous geometries directly.
* **Linear Tetrahedral Meshing:** Generates a 3D structural mesh utilizing 4-node solid tetrahedral elements.
* **Parallel Computing Acceleration:** Dynamically scales computationally intensive individual element stiffness tensor operations across all available local CPU cores via `parfor` parallel loops.
* **Analytical System Inversion:** Reconstructs full system degree-of-freedom (DOF) solutions using high-performance sparse matrix back-substitution solvers.

## 📐 Mathematical Formulation

The linear mechanical behavior of the continuous material domains is governed by classical continuum elasticity equations.

### 1. Kinematic Discretization

For each linear 4-node tetrahedral element, the element volume ($V$) is determined via the determinant of the augmented nodal coordinate matrix:

$$V = \frac{1}{6} \det \begin{bmatrix} 1 & x_1 & y_1 & z_1 \\ 1 & x_2 & y_2 & z_2 \\ 1 & x_3 & y_3 & z_3 \\ 1 & x_4 & y_4 & z_4 \end{bmatrix}$$

The strain-displacement matrix $\mathbf{B}$ maps the local nodal displacements directly into continuous strain tensors utilizing linear shape function spatial derivatives ($\beta_i, \gamma_i, \delta_i$).

### 2. Constitutive Material Matrix

The script maps isotropic material conditions into a $6 \times 6$ elasticity matrix $\mathbf{D}$, defined using Young's modulus ($E$) and Poisson's ratio ($\nu$):

$$\mathbf{D} = \frac{E}{(1+\nu)(1-2\nu)} \begin{bmatrix} 1-\nu & \nu & \nu & 0 & 0 & 0 \\ \nu & 1-\nu & \nu & 0 & 0 & 0 \\ \nu & \nu & 1-\nu & 0 & 0 & 0 \\ 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0 \\ 0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 \\ 0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2} \end{bmatrix}$$

### 3. Numerical Integration & Assembly

Since the strain-displacement relationship is constant across the linear tetrahedron, the local stiffness matrix $\mathbf{k}_e$ simplifies to a direct closed-form evaluation:

$$\mathbf{k}_e = V \cdot (\mathbf{B}^T \mathbf{D} \mathbf{B})$$

The script maps local DOFs to accumulate individual localized entries into a sparse global linear system matrix ($\mathbf{K}_{\text{global}}$):

$$\mathbf{K}_{\text{global}} \mathbf{d} = \mathbf{F}$$

## 💻 Script Configuration Guide

Open the script in MATLAB and calibrate the foundational physical bounds within the `INPUTS` block:

```matlab
%% ================= INPUTS =================
Force = -200;                 % Total compressive force applied at boundary (N)
EE = 2.8e6;                   % Intrinsic Elasticity/Young's modulus (Pa or N/mm^2)
nu = 0.38;                    % Isotropic Poisson's ratio
wu = 3;                       % Axis direction of loading (1=x, 2=y, 3=z)
order_fit = 4;                % Boundary localization tolerance rounding parameter
```

---

# stifness_non_force3.m

**Nonlinear-to-linear equivalent-modulus FEM solver for tetrahedral meshes**, with self-consistent secant extraction and Gibson-Ashby cellular-solids correction.

## Overview

This script solves a static structural finite element problem on an arbitrary tetrahedral mesh (imported from an STL file) under a prescribed applied force. Rather than performing a full incremental nonlinear finite-element solution, it extracts a single **equivalent secant modulus** from a chosen nonlinear constitutive model — calibrated so that a one-shot linear solve reproduces the correct total displacement at the specified applied force — and substitutes it into the standard linear elastic stiffness formulation.

This lets a linear tetrahedral solver reproduce the practical effect of material nonlinearity without the added complexity, convergence sensitivity, and computational cost of a full Newton-Raphson nonlinear solve.

## Supported constitutive models

| # | Model | Description |
|---|-------|-------------|
| 1 | Neo-Hookean | `W = μ/2·(I₁-3)` |
| 2 | Mooney-Rivlin | `W = C10·(I₁-3) + C01·(I₂-3)` |
| 3 | Yeoh | `W = C1·(I₁-3) + C2·(I₁-3)² + C3·(I₁-3)³` |
| 4 | EESM (Arruda-Boyce based) | `W = (1-f)·W_iso + f·W_aniso`, matches the Equivalent Energy Spring Model formulation |
| 5 | Phenomenological | Tangent-modulus power series, `E_tan = E₀·(1+ε+ε/α)` |
| 6 | Polynomial | Direct fit of `Force(λ)`, no strain-energy formulation — for use with `fit()`/`fittype('polyN')` results |

## Algorithm

1. **Strain sweep** builds a `(strain, force)` table by evaluating the chosen model's strain-energy function (or, for the polynomial model, `Force` directly) across a user-specified strain range.
2. **Self-consistent root-find**: solves `F(ε*) = F_applied` exactly via `fzero`, using the swept table only to locate the correct bracket — not a linear-interpolation approximation.
3. **Equivalent modulus extraction**
   - *Linear regime* (small applied strain): the equivalent modulus is taken directly as the independently measured material modulus, `E_t = E`.
   - *Nonlinear regime*: `E_t` is computed as the secant modulus, `E_t = |F_applied| / (A·ε*)` — the physically correct quantity for a one-shot linear solve to reproduce total displacement (as opposed to the local tangent, which represents incremental stiffness rather than displacement at a given load).
4. **Gibson-Ashby correction** (cellular/lattice geometries only): the secant modulus obtained from a lattice specimen's own response is an effective *structural* modulus, not a material property. It is converted to the true constituent material modulus via `E_s = E_lat / (C·(ρ*/ρ_s)^n)` before being used as the mesh's material input — preventing the geometric softening from being counted twice (once implicitly in the secant, once explicitly by the mesh). **This step is skipped for solid specimens** (cube / cube slices), where the secant modulus already is the material modulus.
5. **Linear FEM solve**: assembles the standard isotropic linear-elastic tetrahedral stiffness matrix using `E_t`, solves the reduced system for the applied force, and reconstructs the full displacement field.
6. **Post-processing**: Von Mises stress, displacement contour plots, reaction force.

## Inputs

- STL geometry file, unit system, and scale factor
- Mesh density (`Hmax`) or MATLAB default
- Loading axis and fixed-face selection
- Applied force (sign convention: negative = compression, positive = tension)
- Cross-sectional area `A`
- Constitutive model choice and its parameters (loadable from a saved `.mat` fit file for EESM, or entered manually)
- Strain sweep range and resolution
- Material linear modulus `E` (for the linear-regime check and Gibson-Ashby denominator)
- For lattice geometries: relative density `ρ*/ρ_s`, Gibson-Ashby constant `C`, and exponent `n`

## Outputs

- Self-consistent strain `ε*` and equivalent modulus `E_t`
- Full nodal displacement field
- Von Mises stress distribution
- Maximum displacement and reaction force
- Diagnostic plots: strain-energy sweep, force-strain sweep, geometry, displacement contour, Von Mises contour

## Notes

- The strain-energy sweep and the `fzero` root-find both call the same internal `F_of_eps` function, guaranteeing the swept curve and the self-consistent solution can never disagree with each other.
- For EESM (model 4), the strain-energy formulation matches the corresponding curve-fitting routine (`forceExpression.m`) exactly, so parameters obtained from a `fit()`/`fittype` calibration reproduce the identical force-strain curve when used here.
- This is a **one-shot linear approximation**, not a full nonlinear finite element solution — it reproduces the correct displacement at the specific applied force it was calibrated against, but does not track the full incremental nonlinear path.
