# 3D Linear Static Finite Element Analysis (FEA) for Nodal Displacement

This repository contains a standalone, parallelized MATLAB script for performing 3D linear static structural finite element analysis (FEA). The algorithm imports a customized CAD geometry from an STL file, discretizes the geometry into linear tetrahedral elements, assembles a global stiffness matrix, and evaluates the corresponding 3D displacement vector fields along with boundary reaction forces under localized compressive loading.


## 🛠 Features

* **STL Geometry Ingestion:** Utilizes MATLAB’s Partial Differential Equation (PDE) Toolbox to import and visualize 3D `.stl` continuous geometries directly.
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


%% ================= INPUTS =================
Force = -200;                 % Total compressive force applied at boundary (N)
EE = 2.8e6;                   % Intrinsic Elasticity/Young's modulus (Pa or N/mm^2)
nu = 0.38;                    % Isotropic Poisson's ratio
wu = 3;                       % Axis direction of loading (1=x, 2=y, 3=z)
order_fit = 4;                % Boundary localization tolerance rounding parameter
