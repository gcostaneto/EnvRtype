# Environments in Principal Component Space

Scatter of environments on two PCs, reproducing Fig. 1c of Della Coletta
et al. (2023), where each point is one growth environment labelled by
its code.

## Usage

``` r
env_pca_biplot(
  x,
  pc = c(1, 2),
  label = TRUE,
  group = NULL,
  title = "Environments in PC space",
  verbose = TRUE
)
```

## Arguments

- x:

  an
  [`env_pca`](https://gcostaneto.github.io/EnvRtype/reference/env_pca.md)
  object.

- pc:

  integer vector of length 2. Which PCs to plot. Default `c(1, 2)`.

- label:

  boolean. Draw environment names. Default `TRUE`.

- group:

  factor or named vector. Optional grouping used to colour points, e.g.
  a mega-environment assignment. Names are matched to environment names.

- title:

  character. Plot title.

- verbose:

  boolean. If `TRUE` (default) prints a progress banner.

## Value

A ggplot object when ggplot2 is available (invisibly), otherwise `NULL`
after drawing a base plot.

## See also

[`env_pca`](https://gcostaneto.github.io/EnvRtype/reference/env_pca.md),
[`env_pca_scree`](https://gcostaneto.github.io/EnvRtype/reference/env_pca_scree.md)

## Examples

``` r
# \donttest{
pc <- env_pca(env_indices(maizeWTH, end.day = 130))
#> ---------------------------------------------------------------
#> env_pca -- runs PCA on the environmental index matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ---------------------------------------------------------------
#> env_indices -- builds interval-resolved environmental indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ------------------------------------------------------------
#> env_indices(): 5 environments x 21 factors x 44 intervals = 924 indices
#>   interval width : 3 day(s), days 1..130
#>   discarded      : 106 record(s) outside 1..130
#> ------------------------------------------------------------
#> ------------------------------------------------------------
#> env_pca(): 5 environments x 924 indices
#>   informative PCs: 4 (= n.env - 1)
#>   PC1 54.3%, PC2 21.1%, cumulative 75.4%
#> ------------------------------------------------------------
env_pca_biplot(pc)
#> ---------------------------------------------------------------
#> env_pca_biplot -- plots environments in PC space
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------

env_pca_biplot(pc, pc = c(1, 3))
#> ---------------------------------------------------------------
#> env_pca_biplot -- plots environments in PC space
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------

# }
```
