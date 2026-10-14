# Fit kernel_model() informed by an environment-cluster table

"fixed" adds mega-environment means as fixed effects; "block" zeroes
between-cluster entries of the GxE kernel so genotypes may rank
differently across mega-environments; "both" does each.

## Usage

``` r
kernel_model_clustered(
  y,
  data,
  random,
  env,
  gid,
  clusters,
  use = c("block", "fixed", "both"),
  gxe_pattern = "^KGE_|GE",
  fixed = NULL,
  ...
)
```

## Arguments

- y:

  character. Name of the response column in `data`.

- data:

  data.frame with the response, environment and genotype columns.

- random:

  named list of kernels passed to `kernel_model`.

- env:

  character. Name of the environment column in `data`.

- gid:

  character. Name of the genotype column in `data`.

- clusters:

  an environment-cluster mapping: a `cluster_environments` result, a
  data.frame with columns `env` and `cluster`, or a named vector.

- use:

  character. How clusters enter the model: `"block"` (default),
  `"fixed"` or `"both"`.

- gxe_pattern:

  character. Regular expression selecting the GxE kernels to
  block-diagonalise. Default `"^KGE_|GE"`.

- fixed:

  optional fixed-effects specification passed to `kernel_model`.

- ...:

  further arguments passed to `kernel_model`.

## Value

An object of class `"kernel_model"` (see
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md))
with an added `cluster_info` element: a list recording the
environment-to-cluster `mapping`, the number of clusters `k`, the
`mechanism` used (`"block"`, `"fixed"` or `"both"`), the `gxe_pattern`,
the number of added fixed columns `n_fixed_cols` and the names of the
`blocked_kernels`.

## Details

Block-diagonalisation MODIFIES the GxE kernels, so any cached
decomposition for those kernels is dropped and recomputed. Both dense
and factored kernels are supported. Selected factored kernels are
reconstructed as dense observation-by-observation matrices before
blocking; unselected kernels are left unchanged.
