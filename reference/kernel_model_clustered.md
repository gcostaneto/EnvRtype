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

## Details

Block-diagonalisation MODIFIES the GxE kernels, so any cached
decomposition for those kernels is dropped and recomputed.
