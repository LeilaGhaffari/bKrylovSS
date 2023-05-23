# Shifted Linear Systems:

$$
(\boldsymbol A + \sigma_i \boldsymbol I) \boldsymbol x(\sigma_i) = \boldsymbol b(\sigma_i)
$$

where

$\boldsymbol A \in \mathbb C^{n \times n}$

$i = 1, 2, ..., L$

$\{\sigma_i\}_{i=1}^L \subset \mathbb C$

## Applications:
Tikhonov–Phillips regularization, lattice quantum chromodynamics, rational Krylov subspaces, diﬀuse optical tomography, etc

## Possible Solution:
When $\boldsymbol A$ is *large* and *sparse*, matrix-free iterative methods are of interest, including **Krylov** subspace methods.

### Why Krylov?
The Krylov subspace ($\mathcal{K}_j$) is *invariant* (does not change) under $\boldsymbol A + \sigma_i \boldsymbol I$ (scalar shift of the coeﬃcient matrix).

$$
\mathcal{K}_j(\boldsymbol A + {\sigma_i}_1 \boldsymbol I , \boldsymbol u) = \mathcal{K}_j(\boldsymbol A + {\sigma_i}_2 \boldsymbol I , \widetilde{\boldsymbol u}) \leftarrow Shift \, Invariance
$$

where

$\widetilde{\boldsymbol u} = \beta \boldsymbol u$ ($\beta \neq 0$) $\leftarrow Collinearity \, Requirement$

### Obstacles:

Unrelated RHSs (${\boldsymbol b(\sigma_i)}_{i=1}^L$ not collinear) $\rightarrow$ can't exploit the *shift invariance* property of Krylov subspace

[22] $\rightarrow$ for GMRES-type methods, we cannot simultaneously minimize all residuals while maintaining collinearity (ToDo)

[22] $\rightarrow$ the collinearity requirement causes great difficulty when incorporating shifted system solvers into the subspace recycling framework (?)

### Soodhalter's proposal

Use Krylov subspace:
* shift invariance $\checkmark$
* Collinearity restriction $\times$

But how? Multiple shifted systems is equivalent to a Sylvester equation, therefore the procedure for solving Sylvester equations is exploited [46].

Sylvester operator [17, 46, 55, 53] $\rightarrow$ build bKrylov (block Krylov subspace) $\rightarrow$ shift invariance

[48, 50, 44, 11, 32, 40, 45]

Note!
This paper is not an extension of shifted GMRES [22] or shifted FOM (full orthogonalization method) [54] to bKrylov with multiple RHSs. These extensions already exists [12, 66]. These methods require that the columns of the block residual span the same subspace-hence not applicable to the case of unrelated RHSs.