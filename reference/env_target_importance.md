# Environmental drivers of phenotypic/target clustering

Groups environments by their target metrics (`Y`: phenotypic means,
heritabilities, selection accuracies, etc.) using K-means, then trains a
Random Forest to predict those clusters from the environmental
covariable matrix (`W`). The Random Forest variable importance ranks
which environmental variables best explain the target-based grouping of
environments. The number of clusters can be supplied by the user or
optimised automatically with the average silhouette width.

## Usage

``` r
env_target_importance(
  Y,
  W,
  k = NULL,
  k.max = 10,
  scale.Y = TRUE,
  ntree = 1000,
  importance.type = c("gini", "permutation"),
  nstart = 25,
  seed = NULL,
  save = FALSE,
  dir.path = NULL,
  verbose = TRUE
)
```

## Arguments

- Y:

  numeric vector, matrix or data.frame. Target metrics with environments
  in rows (e.g. adjusted means, \\H^2\\, predictive ability). A named
  vector is also accepted.

- W:

  matrix or data.frame. Environmental covariables with environments in
  rows and variables in columns (e.g. a `W_matrix` output).

- k:

  integer or NULL. Number of clusters. If `NULL` (default) it is chosen
  by maximising the average silhouette width over `2:k.max`.

- k.max:

  integer. Maximum number of clusters tested when `k` is `NULL`. Default
  10.

- scale.Y:

  boolean. If `TRUE` (default) `Y` is centred and scaled before
  clustering.

- ntree:

  integer. Number of trees in the Random Forest. Default 1000.

- importance.type:

  character. `"gini"` (MeanDecreaseGini, default) or `"permutation"`
  (MeanDecreaseAccuracy).

- nstart:

  integer. Number of random starts for K-means. Default 25.

- seed:

  integer or NULL. Optional seed for reproducibility.

- save:

  boolean. If `TRUE`, writes the importance table and cluster assignment
  as `.csv`.

- dir.path:

  character. Output directory when `save = TRUE`. Defaults to the
  working directory.

- verbose:

  boolean. If `TRUE` (default) prints progress messages.

## Value

A list with:

- `clusters`: named integer vector of the cluster of each environment.

- `k`: the number of clusters used.

- `silhouette`: data.frame of average silhouette width per tested `k`
  (`NULL` if `k` was given).

- `kmeans`: the `kmeans` object.

- `randomForest`: the fitted `randomForest` object.

- `importance`: data.frame of environmental variables ranked by
  importance (descending).

## Details

This is a supervised counterpart of `env_cluster`: instead of grouping
environments by their enviromic profile, environments are grouped by
*what the breeder cares about* (yield level, heritability, accuracy) and
the enviromic data are then asked to explain that grouping. Variables at
the top of `importance` are candidate drivers of the observed
environmental stratification and good candidates for
`env_threshold_gxe`.

Permutation importance (`"permutation"`) is generally preferred when
covariables differ in scale or in number of categories, since
MeanDecreaseGini is biased towards high-cardinality predictors.

## References

Costa-Neto, G., da Matta, D., Fernandes, I. K., & Heinemann, A. B.
(2023). Environmental clusters defining breeding zones for tropical
irrigated rice in Brazil. *Agronomy Journal*, 115(5).
[doi:10.1002/agj2.21481](https://doi.org/10.1002/agj2.21481)

## See also

`env_cluster`, `W_matrix`

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
if (requireNamespace("randomForest", quietly = TRUE)) {
  data("maizeYield"); data("maizeWTH")

  ## Environment-level target: mean grain yield per environment
  Y <- tapply(maizeYield$value, maizeYield$env, mean)

  ## Environmental covariables (environments x covariables)
  W <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
                var.id = c("T2M", "T2M_MAX", "PRECTOT"), statistic = "mean")

  ## Which environmental variables best explain the target-based grouping?
  out <- env_target_importance(Y = Y, W = W, seed = 1)
  out$k
  out$clusters
  head(out$importance)

  ## Permutation importance and a fixed number of clusters
  out2 <- env_target_importance(Y = Y, W = W, k = 2,
                                importance.type = "permutation", seed = 1)
  head(out2$importance)

  ## Multi-trait target (e.g. mean and heritability per environment)
  Ymat <- cbind(mean = Y, sd = tapply(maizeYield$value, maizeYield$env, stats::sd))
  out3 <- env_target_importance(Y = Ymat, W = W, seed = 1)
  out3$importance
}
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W
#> ---------------------------------------------------------------
#> env_target_importance -- ranks covariables by importance to a target
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Optimal number of target clusters (silhouette): k = 2
#>   - clustering environments on the target metrics
#>   - fitting random forest of clusters on covariables
#> ---------------------------------------------------------------
#> env_target_importance -- ranks covariables by importance to a target
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - clustering environments on the target metrics
#>   - fitting random forest of clusters on covariables
#> ---------------------------------------------------------------
#> env_target_importance -- ranks covariables by importance to a target
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Optimal number of target clusters (silhouette): k = 2
#>   - clustering environments on the target metrics
#>   - fitting random forest of clusters on covariables
#>       variable importance
#> 1 T2M_MAX_mean  0.8119000
#> 2 PRECTOT_mean  0.6787333
#> 3     T2M_mean  0.6397667
# }
```
