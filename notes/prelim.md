# Preliminaries
1. Krylov subspace methods
1. Techniques for solving Shifted Linear Systems
1. Subspace recycling techniques

## Krylov subspace methods

### Krylov subspace for unshifted system

$\mathbf{Ax} = \mathbf b \xrightarrow[\texttt{ON basis}]{\texttt{Arnoldi}} \mathcal{K}_j(\mathbf A , \mathbf u) = \texttt{span} \{ \mathbf u, \mathbf{Au}, ..., \mathbf A^{j-1} \mathbf u \}$

<small>ON = Orthonormal</small>

### Arnoldi relation:
$\mathbf{A} \mathbf{V}_j = \mathbf{V}_{j+1} \mathbf{\overline{H}}_j$

where

$\mathbf{\overline{H}}_j \in \mathbb C^{(j+1) \times j}$ (upper Hessenberg matrix)

$\mathbf{\overline{V}}_j \in \mathbb C^{n \times j}$ (matrix with orthonormal columns)

$\mathbf x_0 := $  initial approximation to the solution (maybe random)

$\mathbf r_0 = b - \mathbf A \mathbf x_0 := $  initial residual

$\mathbf x_j = \mathbf x_0 + \mathbf t_j$ @ iteration $j$, where $\mathbf t_j \in  \mathcal{K}_j (\mathbf A, \mathbf r_0)$:

$\mathbf t_j :=$ *correction*

$\mathcal{K}_j (\mathbf A, \mathbf r_0) := $ *search space*

### GMRES:

$\mathbf b - \mathbf A (\mathbf x_0 + \mathbf t_j) \bot  \mathbf A \mathcal{K}_j (\mathbf A, \mathbf r_0)$

$\mathbf t_j = \underset{\mathbf t \in \mathcal{K}_j (\mathbf A, \mathbf r_0)}{\texttt{argmin}}  \lVert \mathbf b - \mathbf A (\mathbf x_0 + \mathbf t) \rVert \rightarrow$ this is equivalent to solving a $(j+1) \times j$ minimization problem:

$\mathbf y_j = \underset{\mathbf y \in \mathbb C^j}{\texttt{argmin}}  \left\lVert \overline{H}_j \mathbf y - \left\lVert \mathbf r_0 \right\rVert  \mathbf e_1^{(j+1)} \right\rVert$ $\rightarrow$ $\mathbf x_j = \mathbf x_0 + \mathbf V_j \mathbf y_j$

where $\mathbf e_J^{(i)} :=$ the $J$ th Cartisian basis vector in $\mathbb C^i$

GMRES(m) = `restarted GMRES` : run an $m$ -step cycle of the GMRES method and compute and approximation $\mathbf x_m$

### Full Orthogonalization Method (FOM)

$\mathbf b - \mathbf A (\mathbf x_0 + \mathbf t_j) \bot  \mathbf A \mathcal{K}_j (\mathbf A, \mathbf r_0)$

which is equivalent to solving a $j \times j$ linear system $\mathbf H_j \mathbf y_j = \beta \mathbf e_1^{(j)}$

and where $\mathbf H_j \in \mathbb C^{j \times j}$ is obtained by deleting the last row of $\mathbf{\overline{H}}_j$

FOM(m) = `restarted FOM`

```
Note!
The focus of the paper is the augmentation techniques, Subspace recycling with GMRES (rGMRES or GCRO-DR).

## Recycling GMRES [44]
