# Attach spectral decompositions to a list of dense kernels

Turns each dense kernel (`list(Kernel = <n x n>, Type = "D"|"BD")`) into
a cached spectral decomposition stored on `$decomp`. Call once and reuse
across any number of
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
fits: kernels depend only on genotypes and environments – never on the
phenotype vector – so masking `y` for a cross-validation fold does NOT
invalidate the decomposition.

## Usage

``` r
decompose_kernels(
  K,
  ne = NULL,
  tol = 1e-10,
  scale_kernels = TRUE,
  keep_var = 1,
  verbose = TRUE
)
```

## Arguments

- K:

  list of kernels.

- ne:

  per-environment observation counts; defaults to `attr(K, "ne")`.

- tol:

  eigenvalue cutoff; must match `kernel_model(tol=)`.

- scale_kernels:

  must match `kernel_model(scale_kernels=)`.

- keep_var:

  numeric in (0, 1\]. Retain only the leading eigenvalues carrying this
  proportion of each kernel's trace. 1 (default) keeps everything. A
  named vector or list keyed by kernel name / regex applies per kernel,
  e.g. `c(KGE_ = 0.95)`.

- verbose:

  report rank, timing and route per kernel.

## Value

`K` with `$decomp` on each element and a `"kd_meta"` attribute.

## Details

Factored kernels from
[`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md)
arrive pre-decomposed and do not need this function.

## See also

[`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md),
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md),
[`undecompose_kernels`](https://gcostaneto.github.io/EnvRtype/reference/undecompose_kernels.md)
