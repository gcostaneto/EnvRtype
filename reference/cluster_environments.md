# Group environments into mega-environments (clusters)

Delineates mega-environments by clustering the **genetic correlation**
among environments, i.e. how similarly genotypes rank from one
environment to another. Environments in which genotypes respond alike
are placed in the same cluster, which is the operational definition of a
mega-environment in a target population of environments (TPE).

The input is either a fitted
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
(carrying the environment genetic-correlation matrix in
`$genetic_cor_env`, produced with `keep_effects = TRUE`) or any square
environment correlation matrix you already have (for example from
[`env_cor`](https://gcostaneto.github.io/EnvRtype/reference/env_cor.md)
or an external GxE analysis).

## Usage

``` r
cluster_environments(
  object,
  k = NULL,
  method = c("hclust", "pam", "kmeans"),
  distance = c("one_minus", "sqrt"),
  hclust_linkage = "ward.D2",
  k_range = 2:6,
  seed = NULL
)
```

## Arguments

- object:

  a
  [`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
  fit carrying `$genetic_cor_env` (fit with `keep_effects = TRUE`), or a
  square environment-by-environment correlation matrix.

- k:

  integer or NULL. Number of clusters. When `NULL` (default) it is
  chosen by maximum average silhouette width over `k_range`.

- method:

  character. Clustering algorithm: `"hclust"` (default), `"pam"` or
  `"kmeans"`.

- distance:

  character. Distance built from the genetic correlation: `"one_minus"`
  (\\1 - r_g\\, default) or `"sqrt"` (\\\sqrt{1 - r_g}\\).

- hclust_linkage:

  character. Linkage passed to
  [`hclust`](https://rdrr.io/r/stats/hclust.html) when
  `method = "hclust"`. Default `"ward.D2"`.

- k_range:

  integer vector. Candidate numbers of clusters tested when `k` is
  `NULL`. Default `2:6`.

- seed:

  integer or NULL. Optional seed for reproducibility of the randomised
  starts in `"pam"` and `"kmeans"`.

## Value

A list with:

- clusters:

  named integer vector giving the cluster of each environment.

- k:

  the number of clusters used (selected or supplied).

- method:

  the clustering algorithm used.

- distance:

  the symmetric distance matrix that was clustered.

- silhouette:

  per-environment silhouette widths plus their `average` (higher is a
  tighter, better-separated solution).

- medoids:

  representative environment of each cluster (`"pam"`/`"kmeans"`; `NULL`
  for `"hclust"`).

- hclust:

  the [`hclust`](https://rdrr.io/r/stats/hclust.html) tree for
  `method = "hclust"`, otherwise `NULL`.

- table:

  a data.frame (`env`, `cluster`, `sil_width`) ordered by cluster and
  decreasing silhouette width.

## Details

The genetic correlation \\r_g\\ is turned into a distance and the chosen
algorithm is run on that distance:

- `distance = "one_minus"` uses \\d = 1 - r_g\\ (default).

- `distance = "sqrt"` uses \\d = \sqrt{1 - r_g}\\, which spreads out
  highly correlated environments.

Correlations are clamped to \\\[-1, 1\]\\ and the distance matrix is
symmetrised before clustering. Three algorithms are available through
`method`: agglomerative hierarchical clustering (`"hclust"`, controlled
by `hclust_linkage`), partitioning around medoids (`"pam"`), and
`"kmeans"` run on a classical MDS embedding of the distance matrix.

When `k` is `NULL` the number of clusters is selected automatically by
maximising the **average silhouette width** over the candidate values in
`k_range`; supplying `k` skips this search. At least three environments
are required.

## See also

[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
(fit with `keep_effects = TRUE` to obtain `$genetic_cor_env`),
[`kernel_model_clustered`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model_clustered.md)
to fit a model informed by the clusters,
[`env_cluster`](https://gcostaneto.github.io/EnvRtype/reference/env_cluster.md)
for clustering from an environmental covariable matrix, and
[`env_cor`](https://gcostaneto.github.io/EnvRtype/reference/env_cor.md).

## Examples

``` r
## ------------------------------------------------------------------
## 1. Cluster a genetic-correlation matrix directly (no model needed)
## ------------------------------------------------------------------
## Six environments forming two obvious groups
R <- matrix(0.2, 6, 6)
R[1:3, 1:3] <- 0.9
R[4:6, 4:6] <- 0.85
diag(R) <- 1
dimnames(R) <- list(paste0("E", 1:6), paste0("E", 1:6))

## Automatic number of clusters via silhouette width
cl <- cluster_environments(R)
cl$k
#> [1] 4
cl$table
#>   env cluster sil_width
#> 1  E1       1     0.875
#> 2  E2       1     0.875
#> 3  E3       1     0.875
#> 4  E4       2     1.000
#> 5  E5       3     1.000
#> 6  E6       4     1.000
cl$silhouette["average"]
#> average 
#>  0.9375 

## ------------------------------------------------------------------
## 2. Force a specific number of clusters
## ------------------------------------------------------------------
cluster_environments(R, k = 2)$clusters
#> E1 E2 E3 E4 E5 E6 
#>  1  1  1  2  2  2 

## ------------------------------------------------------------------
## 3. Try different algorithms and distances
## ------------------------------------------------------------------
cluster_environments(R, k = 2, method = "pam", seed = 1)$medoids
#> [1] "E1" "E4"
cluster_environments(R, method = "kmeans", k_range = 2:4, seed = 1)$clusters
#> E1 E2 E3 E4 E5 E6 
#>  3  3  3  1  4  2 
cluster_environments(R, distance = "sqrt",
                     hclust_linkage = "complete")$clusters
#> E1 E2 E3 E4 E5 E6 
#>  1  1  1  2  3  4 

if (FALSE) { # \dontrun{
## ------------------------------------------------------------------
## 4. Full workflow from phenotypes, kernels and a fitted model
## ------------------------------------------------------------------
data("maizeYield"); data("maizeG"); data("maizeWTH")

ECs <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
                var.id = c("FRUE", "PETP", "SRAD", "T2M_MAX"))
K <- get_kernel(K_G = list(G = maizeG),
                K_E = list(W = env_kernel(env.data = ECs)[[2]]),
                data = maizeYield, model = "RNMM")

## keep_effects = TRUE is required for $genetic_cor_env
fit <- kernel_model(y = "value", env = "env", gid = "gid",
                    data = maizeYield, random = K,
                    keep_effects = TRUE)

megaenv <- cluster_environments(fit)
megaenv$table

## 5. Feed the clusters into a cluster-aware model
fit_cl <- kernel_model_clustered(y = "value", data = maizeYield,
                                 random = K, env = "env", gid = "gid",
                                 clusters = megaenv, use = "block")
} # }
```
