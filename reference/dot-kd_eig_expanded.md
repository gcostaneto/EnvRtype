# Exact spectral decomposition of K = M\[idx, idx\] without forming n x n eigen.

M = V D V' =\> K = (V\[idx,\]) D (V\[idx,\])'. Let A = V\[idx,\] sqrt(D)
so K = A A'. The nonzero eigenpairs of A A' follow from the small Gram
matrix A'A = W L W': eigenvalues L, eigenvectors A W L^-1/2. Exact (not
Nystrom): verified to 1e-14 against full eigen().

## Usage

``` r
.kd_eig_expanded(M, idx, tol = 1e-10)
```
