# Build the Matrix of Environmental Covariables for Reaction Norm

Estimates a weather covariable realized matrix **W** with dimensions \\q
\times k\\, for \\q\\ environments and \\k\\ covariables. This matrix is
the basic input of the reaction-norm and enviromic-assisted models in
EnvRtype (`env_kernel`, `env_cluster`, `env_target_importance`,
`get_kernel`).

## Usage

``` r
W_matrix(
  env.data,
  is.processed = FALSE,
  id.names = NULL,
  env.id = NULL,
  var.id = NULL,
  probs = NULL,
  by.interval = NULL,
  time.window = NULL,
  names.window = NULL,
  center = TRUE,
  scale = TRUE,
  sd.tol = 10,
  statistic = NULL,
  tol = 0.001,
  QC = FALSE,
  impute = c("none", "mean", "drop"),
  verbose = TRUE,
  copula = NULL,
  copula.args = list()
)
```

## Arguments

- env.data:

  data.frame of environmental variables from `get_weather`, or an
  already summarised object when `is.processed = TRUE`.

- is.processed:

  boolean. Indicates whether the data.frame was previously processed
  with `summaryWTH` and already contains means, medians, etc.

- id.names:

  character. Columns used as id for the environmental variables.

- env.id:

  character. Column used as id for environments.

- var.id:

  vector (character). Which variables will be used in the analysis.

- probs:

  vector (numeric). Probability quantiles in \\\[0,1\]\\, used when
  `statistic = 'quantile'`. If `NULL`, `c(0.25, 0.50, 0.75)`.

- by.interval:

  boolean. Indicates if temporal intervals must be computed inside each
  environment. Default `FALSE`.

- time.window:

  vector (numeric). If `by.interval = TRUE`, the temporal breaks.

- names.window:

  vector (character). If `by.interval = TRUE`, the interval names.

- center:

  boolean. Indicates whether the matrix should be centred. Default
  `TRUE`.

- scale:

  boolean. If `TRUE`, variables assume \\x \sim N(0,1)\\. Default
  `TRUE`.

- sd.tol:

  numeric. Maximum standard deviation tolerated in quality control.
  Default 10.

- statistic:

  vector (character). Statistic to be computed,
  `c('all','sum','mean','quantile')`. Default `'mean'`.

- tol:

  numeric. Numerical tolerance; variables whose standard deviation is at
  or below `tol` are treated as near-constant and flagged for removal.
  Default 1E-3.

- QC:

  boolean. Indicates whether Quality Control is applied. QC removes
  variables with `sd(x) > sd.tol` and near-constant variables.

- impute:

  character. How to handle missing values before scaling: `"none"`
  (default, keep NAs), `"mean"` (column mean) or `"drop"` (remove
  columns with any NA).

- verbose:

  boolean. If `TRUE` (default) prints quality-control messages.

- copula:

  character or NULL. If not `NULL`, the covariable matrix is replaced by
  copula-based indices computed by `env_copula`. One of `"none"`
  (default, same as `NULL`), `"pobs"`, `"joint"`, `"survival"`,
  `"kendall"` or `"all"`. See `env_copula` for the meaning of each
  index.

- copula.args:

  list. Extra arguments passed to `env_copula`, e.g.
  `list(groups = list(heat = c("T2M_mean","T2M_MAX_mean")), ties.method = "average")`.

## Value

An environmental covariable realized matrix with dimensions \\q \times
k\\. The centring and scaling values, plus the list of removed markers,
are attached as attributes (`"scaled:center"`, `"scaled:scale"`,
`"removed"`) so that new environments can be projected onto the same
space.

## Details

Quality control follows Morais Junior et al. (2018): covariables whose
standard deviation across environments exceeds `sd.tol` are discarded,
as are near-constant covariables (\\sd \le tol\\) which carry no
information about environmental differences.

Unlike earlier versions, the numerical tolerance is *not* added to the
data before scaling (which silently shifted unscaled outputs); it is
used only to detect near-constant columns.

## References

Costa-Neto G., Galli G., Carvalho H.F., Crossa J., Fritsche-Neto R.
(2021). EnvRtype: a software to interplay enviromics and quantitative
genomics in agriculture. *G3* 11(4), jkab040.

## See also

`summaryWTH`, `env_kernel`, `env_typing`,
[`T_matrix`](https://gcostaneto.github.io/EnvRtype/reference/T_matrix.md)

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
data("maizeWTH")
env.data <- maizeWTH[maizeWTH$daysFromStart < 100, ]

## Mean-centred and scaled matrix (default statistic = 'mean')
W <- W_matrix(env.data = env.data)
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W
dim(W)
#> [1]  5 21

## Adding time windows (one block of covariables per development stage)
W <- W_matrix(env.data = env.data, by.interval = TRUE,
              time.window = c(0, 14, 35, 60, 90, 120))
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W

## Selecting the statistic to be used
W <- W_matrix(env.data = env.data, by.interval = TRUE, statistic = 'quantile',
              time.window = c(0, 14, 35, 60, 90, 120))
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W

## With Quality Control based on the maximum sd tolerated
W <- W_matrix(env.data = env.data, QC = TRUE, sd.tol = 3)
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W
#> ------------------------------------------------
#> Quality Control based on sd.tol = 3
#> Removed variables: 2 from 21
#> ALLSKY_TOA_SW_DWN_mean
#> RH2M_mean
#> ------------------------------------------------
attr(W, "removed")
#> [1] "ALLSKY_TOA_SW_DWN_mean" "RH2M_mean"             

## Creating W for specific variables
W <- W_matrix(env.data = env.data, var.id = c('T2M_MAX', 'T2M_MIN', 'T2M'))
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W

## Combining with summaryWTH by using is.processed = TRUE
data <- summaryWTH(env.data, env.id = 'env', statistic = 'quantile')
#> ---------------------------------------------------------------
#> summaryWTH -- summarises weather by environment and interval
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
W <- W_matrix(env.data = data, is.processed = TRUE)
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - centring, scaling and quality-controlling W
#> Error in W_matrix(env.data = data, is.processed = TRUE): object 'W' not found
# }
```
