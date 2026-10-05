# Exact orthonormal decomposition of K = Kb\[idx, idx\]

Never forms the `n x n` matrix. With `Kb = V D V'}, \code{K = A A'` for
`A = V[idx,] sqrt(D)`; the nonzero eigenpairs of
`A A'} follow from the small Gram matrix \code{A'A = W L W'} as eigenvalues \code{L} and eigenvectors \code{A W L^-1/2}. } \details{ This routes through the Gram matrix rather than returning \code{V[idx,]} directly. \code{V[idx,]} is \emph{not} orthonormal when \code{idx} repeats -- its columns have norm \code{sqrt(multiplicity)} -- and although \code{K = U S U'`
still reconstructs exactly,
[`kernel_model()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
uses `U` and `deltav = 1/s` separately, so the error does not cancel.

## Usage

``` r
.fk_expand(Kb, idx, tol = 1e-10, keep_var = 1, scale = TRUE)
```

## Arguments

- Kb:

  small symmetric matrix (`p x p`).

- idx:

  integer vector mapping observations to rows of `Kb`.

- tol:

  relative eigenvalue cutoff.

- keep_var:

  proportion of the trace to retain.

- scale:

  normalise to `mean(diag(K)) = 1`.

## Value

list with `s`, `r`, `n` and `U()`.
