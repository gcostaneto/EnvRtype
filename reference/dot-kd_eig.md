# Eigendecompose a symmetric matrix, dropping negligible eigenvalues. eigen(symmetric = TRUE) returns DECREASING eigenvalues, so the retained set is a leading prefix and we slice rather than which().

Eigendecompose a symmetric matrix, dropping negligible eigenvalues.
eigen(symmetric = TRUE) returns DECREASING eigenvalues, so the retained
set is a leading prefix and we slice rather than which().

## Usage

``` r
.kd_eig(K, tol = 1e-10)
```
