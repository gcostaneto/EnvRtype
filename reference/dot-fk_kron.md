# Exact decomposition of a Hadamard (GxE) kernel via its Kronecker form

The elementwise product of two expansions is the expansion of a
Kronecker product: with `j = (ie-1)*ng + ig`, `K = (Ke %x% Kg)[j, j]`.

## Usage

``` r
.fk_kron(Kg, ig, Ke, ie, tol = 1e-10, keep_var = 1, scale = TRUE)
```

## Arguments

- Kg, Ke:

  small genomic and environmental matrices.

- ig, ie:

  observation index vectors.

- tol, keep_var, scale:

  as in `.fk_expand`.

## Value

list with `s`, `r`, `n` and `U()`.

## Details

On a complete balanced grid the eigenpairs are closed-form – eigenvalues
`outer(dg, de)` and eigenvectors the Khatri-Rao product of rows – so no
large eigendecomposition is needed and truncation is applied *before*
the eigenvectors are materialised. On unbalanced designs that
construction is not orthonormal, so the Gram route is used instead.
