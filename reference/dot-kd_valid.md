# Are cached decompositions valid for the requested settings?

A mismatch in `tol` or `scale_kernels` silently changes the model, so
[`kernel_model()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
re-decomposes rather than trust a stale cache.

## Usage

``` r
.kd_valid(K, tol, scale_kernels, n)
```

## Arguments

- K:

  a kernel list.

- tol, scale_kernels, n:

  settings to validate against.

## Value

logical.
