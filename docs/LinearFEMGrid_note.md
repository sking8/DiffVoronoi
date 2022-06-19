## Linear FEM Grid

The basic mathematical formula for linear FEM Grid is 
$$
\begin{equation}
    \bm K \bm u=\bm f
\end{equation}
$$
where $K$ is the global stiffness matrix consisting of element stiffness matricies. This equation is the main focus of `LinearFEMGrid.h`.

The $\bm K^{(e)}$ is calculated in `LinearFEMFunc.h`. Each $\bm K^{(e)}$ is uniquely determined by the material parameters `(youngs modulus, poisson's ratio)`.
$$
\begin{equation}
    \bm K^{(e)}=\int_{V^{(e)}} \bm{B^T}\bm E\bm B dV^{(e)}
\end{equation}
$$

Here $\bm E$ is the stress-strain matrix and $B$ is the strain-displacement matrix. After discretization using Gauss quadrature rules (assuming the 3D case), the formula becomes:

$$
\begin{equation}
    \bm K^{(e)} = \sum_{i=1}^{p_1}\sum_{j=1}^{p_2}\sum_{k=1}^{p_3}w_i w_j w_k \bm{B^T_{ijk}}\bm E\bm B_{ijk}J_{ijk}
\end{equation}
$$
$J$ in code is `dNde`, meaning the derivatives of shape function w.r.t. Cartetisan coordinates.

Other variables in code:

G is the shear modulus, a.k.a. modulus of rigidity. It has physical dimension: stress=force/area (e.g. MPa).

e is the effective modulus modified by Poisson's ratio