# Coerce a Simulated Envirome to a Plain Matrix

Strips the `sim_W` class and its diagnostic attributes, returning the
bare \\q \times n\_{var}\\ covariable matrix for use with
[`env_kernel`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md),
[`W_matrix`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md)
consumers, or base matrix ops.

## Usage

``` r
# S3 method for class 'sim_W'
as.matrix(x, ...)
```

## Arguments

- x:

  a `sim_W` object.

- ...:

  ignored.

## Value

A numeric matrix (environments in rows, covariables in columns) with no
`sim_W` class or attributes.

## See also

[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md)
