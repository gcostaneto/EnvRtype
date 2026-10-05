# Envirotype clustering into mega-environments

Clusters environments into homogeneous groups ("mega-environments") from
an environmental covariable matrix (`W`) using K-means. The number of
groups can be set by the user or selected automatically by maximising
the average silhouette width. An environmental relationship (kinship)
matrix is also returned for convenience.

## Usage

``` r
env_cluster(
  W,
  k = NULL,
  k.max = 10,
  scale = TRUE,
  nstart = 25,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- W:

  matrix or data.frame. Environmental covariables (environments in rows,
  variables in columns), typically a `W_matrix` output. An envirotype
  frequency matrix from
  [`T_matrix`](https://gcostaneto.github.io/EnvRtype/reference/T_matrix.md)
  works identically.

- k:

  integer or NULL. Number of mega-environments. If `NULL` (default) it
  is chosen by silhouette.

- k.max:

  integer. Maximum number of clusters tested when `k` is `NULL`. Default
  10.

- scale:

  boolean. If `TRUE` (default) `W` is centred and scaled before
  clustering.

- nstart:

  integer. Number of random starts for K-means. Default 25.

- seed:

  integer or NULL. Optional seed for reproducibility.

- verbose:

  boolean. If `TRUE` (default) prints progress messages.

## Value

A list with `clusters` (named integer vector), `k`, `silhouette`
(data.frame of average silhouette per tested `k`, or `NULL` when `k` was
supplied), `kmeans` (the fitted object) and `env.kinship` (an
environment x environment correlation matrix).

## Details

Mega-environments are repeatable groups of locations that rank genotypes
similarly. Here they are approximated by grouping environments that
share a similar enviromic profile, which is useful to stratify a Target
Population of Environments (TPE) before multi-environment prediction.
The silhouette width \\s(i)\\ contrasts within- and between-cluster
distances; the `k` maximising its average is retained.

## See also

`W_matrix`, `env_kernel`, `env_target_importance`,
[`T_matrix`](https://gcostaneto.github.io/EnvRtype/reference/T_matrix.md)

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
data("maizeWTH")

W <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
              var.id = c("T2M", "T2M_MAX", "PRECTOT"), statistic = "mean")
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W

## Automatic number of mega-environments (chosen by silhouette)
mega <- env_cluster(W, seed = 1)
#> ---------------------------------------------------------------
#> env_cluster -- groups environments into mega-environments
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - selecting the number of clusters by silhouette
#> Optimal number of mega-environments (silhouette): k = 3
#>   - running k-means clustering
mega$k
#> [1] 3
mega$clusters
#> NM SO PM IP SE 
#>  3  2  1  3  1 
mega$silhouette
#>   k avg.silhouette
#> 1 2      0.4278887
#> 2 3      0.5017074
#> 3 4      0.2443995

## Force a fixed number of groups
mega2 <- env_cluster(W, k = 2, seed = 1)
#> ---------------------------------------------------------------
#> env_cluster -- groups environments into mega-environments
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - running k-means clustering
mega2$clusters
#> NM SO PM IP SE 
#>  1  2  1  1  1 

## Environmental relationship (kinship) matrix
round(mega$env.kinship, 2)
#>       NM    SO    PM    IP    SE
#> NM  1.00 -0.73  0.35  0.72  0.56
#> SO -0.73  1.00 -0.90 -1.00 -0.98
#> PM  0.35 -0.90  1.00  0.90  0.97
#> IP  0.72 -1.00  0.90  1.00  0.98
#> SE  0.56 -0.98  0.97  0.98  1.00

## Use the grouping downstream (e.g. as a fixed effect or to stratify CV)
table(mega$clusters)
#> 
#> 1 2 3 
#> 2 1 2 
# }
```
