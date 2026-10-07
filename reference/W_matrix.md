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
  impute = c("none", "mean", "median", "knn", "drop"),
  knn = 5L,
  max.cor = NULL,
  group = FALSE,
  cor.report = NULL,
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
  (default, keep NAs), `"mean"` (column mean), `"median"` (column
  median), `"knn"` (k-nearest-environment imputation, see `knn`) or
  `"drop"` (remove columns with any NA).

- knn:

  integer. Number of nearest environments used when `impute = "knn"`.
  Each missing cell is filled with the mean of that covariable over the
  `knn` most similar environments (Euclidean distance on the
  standardised, commonly observed covariables). Default 5.

- max.cor:

  numeric in \\(0,1\]\\ or `NULL`. If not `NULL`, a greedy collinearity
  filter (same rule as `caret::findCorrelation`) drops covariables so
  that no pair of retained covariables has absolute Pearson correlation
  above `max.cor`. When two covariables are too correlated the one with
  the larger mean absolute correlation to the rest is removed. Unlike
  the `sd`-based rules this is an explicit opt-in and the collinear
  covariables are always dropped (independently of `QC`). Default `NULL`
  (no collinearity filtering).

- group:

  boolean. Only used when `max.cor` is not `NULL`. If `TRUE`, collinear
  covariables are *grouped* instead of dropped: covariables are
  clustered by average-linkage hierarchical clustering on \\1 - \|r\|\\,
  the tree is cut at height \\1 - \\`max.cor`, and each correlated block
  is replaced by its mean (a single composite covariable). The
  block-to-member mapping is returned in the `"groups"` attribute.
  Default `FALSE` (collinear covariables are dropped).

- cor.report:

  `NULL`, `TRUE`, or a character path. If not `NULL`, a per-covariable
  collinearity diagnosis (its largest absolute correlation, the partner
  responsible, the final status, and the block it belongs to) is written
  as a CSV. Pass `TRUE` to write `"W_matrix_collinearity.csv"` in the
  current working directory, a directory to write that file there, or a
  full file path. Default `NULL` (no file written). The same table is
  always attached as the `"collinearity"` attribute.

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
space. When quality control or collinearity filtering drop covariables,
the reason for each removal is reported in the `"removed.reason"`
attribute (a named character vector with values `"too.variable"`,
`"near.constant"` or `"collinear"`). The full collinearity diagnosis is
attached as the `"collinearity"` data.frame attribute, and when
`group = TRUE` the block-to-member mapping is attached as `"groups"`.

## Details

Quality control follows Morais Junior et al. (2018): covariables whose
standard deviation across environments exceeds `sd.tol` are discarded,
as are near-constant covariables (\\sd \le tol\\) which carry no
information about environmental differences.

Three further, independent cleaning steps are available. Missing values
can be imputed column-wise (`impute = "mean"`/`"median"`) or from the
most similar environments (`impute = "knn"`) before any statistic is
computed, which avoids silently shrinking the covariable set the way
`impute = "drop"` does. Redundant covariables can be pruned with
`max.cor`: enviromic matrices frequently contain near-duplicated columns
(e.g. several temperature summaries), which inflate Euclidean geometry
and relatedness kernels; the greedy filter keeps one representative per
correlated block.

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
#> Quality Control (sd.tol = 3)
#> Removed variables: 2 from 21
#>   too variable (sd > sd.tol): ALLSKY_TOA_SW_DWN_mean, RH2M_mean
#> ------------------------------------------------
attr(W, "removed")
#> [1] "ALLSKY_TOA_SW_DWN_mean" "RH2M_mean"             

## Dropping redundant (collinear) covariables and imputing missing cells
W <- W_matrix(env.data = env.data, impute = "knn", max.cor = 0.95)
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W
#> ------------------------------------------------
#> Quality Control (sd.tol = 10, max.cor = 0.95)
#> Removed variables: 9 from 21
#>   collinear (|r| > max.cor): GDD_mean, T2M_mean, FRUE_mean, RTA_mean, ETP_mean, T2M_MIN_mean, SPV_mean, T2MDEW_mean, ALLSKY_SFC_SW_DWN_mean
#> ------------------------------------------------
attr(W, "removed.reason")
#>               GDD_mean               T2M_mean              FRUE_mean 
#>            "collinear"            "collinear"            "collinear" 
#>               RTA_mean               ETP_mean           T2M_MIN_mean 
#>            "collinear"            "collinear"            "collinear" 
#>               SPV_mean            T2MDEW_mean ALLSKY_SFC_SW_DWN_mean 
#>            "collinear"            "collinear"            "collinear" 

## Grouping collinear covariables into composite blocks + CSV diagnosis
W <- W_matrix(env.data = env.data, max.cor = 0.9, group = TRUE,
              cor.report = tempdir())
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W
#> ------------------------------------------------
#> Collinearity grouping (max.cor = 0.9)
#> 21 covariables -> 8 blocks
#> ------------------------------------------------
#> Collinearity diagnosis written to: /tmp/Rtmp2OrnoG/W_matrix_collinearity.csv
attr(W, "groups")
#>       ALLSKY_SFC_LW_DWN_mean       ALLSKY_SFC_SW_DWN_mean 
#> "ALLSKY_SFC_LW_DWN_mean(+1)" "ALLSKY_TOA_SW_DWN_mean(+2)" 
#>       ALLSKY_TOA_SW_DWN_mean                     ETP_mean 
#> "ALLSKY_TOA_SW_DWN_mean(+2)"               "ETP_mean(+1)" 
#>                    FRUE_mean                     GDD_mean 
#>               "T2M_mean(+6)"               "T2M_mean(+6)" 
#>                       N_mean                    PETP_mean 
#>                 "N_mean(+1)"              "PETP_mean(+1)" 
#>                 PRECTOT_mean                    RH2M_mean 
#>              "PETP_mean(+1)"              "RH2M_mean(+1)" 
#>                     RTA_mean                     SPV_mean 
#>                 "N_mean(+1)"               "T2M_mean(+6)" 
#>                    SRAD_mean                  T2MDEW_mean 
#>               "ETP_mean(+1)" "ALLSKY_SFC_LW_DWN_mean(+1)" 
#>                 T2M_MAX_mean                 T2M_MIN_mean 
#>               "T2M_mean(+6)"               "T2M_mean(+6)" 
#>               T2M_RANGE_mean                     T2M_mean 
#> "ALLSKY_TOA_SW_DWN_mean(+2)"               "T2M_mean(+6)" 
#>                     VPD_mean                    WS2M_mean 
#>              "RH2M_mean(+1)"               "T2M_mean(+6)" 
#>                       n_mean 
#>                     "n_mean" 

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
