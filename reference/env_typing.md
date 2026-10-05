# Environmental Typologies based on Cardinal, Quantilic Limits or Data-Driven Mining

Returns environmental typologies (*envirotypes*) that can be used as
envirotype markers. Typologies are given by cardinals (discrete
intervals for each variable), by empirical quantiles, or – when
`envirotype_mining = TRUE` – learned from the data: for each variable a
K-means clustering is run over a grid of candidate `k` values and the
optimal number of envirotypes is selected by the Calinski-Harabasz
(variance ratio) criterion.

## Usage

``` r
env_typing(
  env.data,
  var.id,
  env.id,
  cardinals = NULL,
  days.id = NULL,
  time.window = NULL,
  names.window = NULL,
  quantiles = NULL,
  id.names = NULL,
  by.interval = FALSE,
  scale = FALSE,
  format = NULL,
  ratio = FALSE,
  joint = FALSE,
  joint.by.interval = TRUE,
  combine.features = FALSE,
  combine.max.order = 2L,
  envirotype_mining = FALSE,
  k.range = 2:10,
  nstart = 25,
  iter.max = 100,
  seed = 1234,
  plot_envirotypes = FALSE,
  min.window = 2L,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame of environmental variables from `get_weather`.

- var.id:

  character. Which variables will be used in the analysis.

- env.id:

  character. Column used as id for environments.

- cardinals:

  list (numeric). Cardinal thresholds per variable. If `NULL`, see
  `quantiles`. Ignored when `envirotype_mining = TRUE`.

- days.id:

  character. Name of the column indicating the days from start.

- time.window:

  vector (numeric). If `by.interval = TRUE`, the temporal breaks.

- names.window:

  vector (character). If `by.interval = TRUE`, the interval names. Note
  that [`cut()`](https://rdrr.io/r/base/cut.html) is applied with an
  open-ended final break, so `time.window = c(0,40,80,120)` defines FOUR
  intervals; supply either 4 names or 3, in which case the trailing one
  is auto-named.

- quantiles:

  vector (numeric). Probability quantiles used when `cardinals` is
  `NULL`. Default `c(0.01, .25, .50, .75, .99)`. Ignored when mining.

- id.names:

  vector (character). Columns used as id for the environmental
  variables.

- by.interval:

  boolean. Compute temporal intervals inside each environment. Default
  `FALSE`.

- scale:

  boolean. If `TRUE`, variables assume \\x \sim N(0,1)\\. Default
  `FALSE`.

- format:

  character. Output shape, `'long'` (default) or `'wide'`.

- ratio:

  boolean. If `TRUE`, `Freq` is rescaled to a relative frequency in
  \\\[0,1\]\\ within each environment x interval x variable group.
  Default `FALSE`.

- joint:

  boolean. Only used when `envirotype_mining = TRUE`. If `FALSE`
  (default) one univariate K-means is run per variable. If `TRUE`, the
  variables in `var.id` are COMBINED: the vector of all variables
  observed at the same time point is clustered jointly, so each
  envirotype is a recurring multi-variable state of the environment
  (analogous to a haplotype across loci).

- joint.by.interval:

  boolean. Only used when `joint = TRUE`. If `TRUE` (default) a separate
  K-means is fitted inside each time interval; if `FALSE` a single
  K-means is fitted over all intervals pooled.

- combine.features:

  boolean. Only valid when `envirotype_mining = TRUE` and
  `joint = FALSE`. If `TRUE`, envirotypes are mined for every unique
  COMBINATION of the variables in `var.id`: each variable on its own,
  then every pair, every triple, and so on. A combination of two or more
  variables is clustered jointly (as in `joint = TRUE`); a
  single-variable combination reduces to the univariate case. Each
  combination receives a compact generic id (`FC001`, `FC002`, ...) and
  the returned `feature_combinations` table maps every id back to its
  variables and time intervals. Default `FALSE`.

- combine.max.order:

  integer or `NULL`. Largest combination size to enumerate when
  `combine.features = TRUE`. Because the number of combinations grows as
  a power set (\\2^n - 1\\), this defaults to `2` (single variables and
  pairs only), which stays tractable for dozens of variables. Increase
  it to include triples and higher (e.g. `3`), or set `NULL` to
  enumerate all `length(var.id)` orders. A guard stops the run if the
  requested number of combinations is impractically large.

- envirotype_mining:

  boolean. If `TRUE`, envirotype classes are mined from the data by
  K-means with `k` chosen by the Calinski-Harabasz index. Default
  `FALSE`.

- k.range:

  integer vector. Candidate numbers of envirotypes tested when mining.
  Default `2:10`.

- nstart:

  numeric. Random starts passed to
  [`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html). Default 25.

- iter.max:

  numeric. Maximum iterations passed to
  [`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html). Default 100.

- seed:

  numeric. Random seed for reproducible K-means. Default 1234; `NULL`
  skips seeding.

- plot_envirotypes:

  boolean. If `TRUE`, draws a barplot of the frequency (or ratio) of
  each envirotype per environment, faceted by variable (univariate) or
  by time interval (joint). Uses ggplot2 when available, base graphics
  otherwise. Default `FALSE`.

- min.window:

  integer. Minimum number of observations required in the open-ended
  trailing interval. Default 2; use 0 or `NULL` to always keep it.

- verbose:

  boolean. If `TRUE` (default) prints progress messages.

## Value

If `envirotype_mining = FALSE` (default) and `plot_envirotypes = FALSE`:
a data.frame (or matrix, if `format = 'wide'`) of envirotype
frequencies.

Otherwise a `list` with:

- `typologies`:

  data.frame (long) or matrix (wide) of envirotype frequencies.

- `envirotype_description`:

  data.frame describing the numeric range captured by each envirotype:
  variable, label, cluster index, `n`, `min`, `max`, `mean`, `sd`,
  `median` and the `lower`/`upper` cut boundaries.

- `mining_summary`:

  data.frame with the optimal `k`, the CH index at the optimum and the
  full CH profile over `k.range`.

- `feature_combinations`:

  only when `combine.features = TRUE`: a data.frame mapping each generic
  combination id to its number of features (`n_features`), the
  `features` themselves and the time `intervals`.

- `plot`:

  the plot object, or `NULL`.

## Details

The Calinski-Harabasz index for a partition into \\k\\ groups of \\n\\
observations is \$\$CH(k) = \frac{BGSS/(k-1)}{WGSS/(n-k)}\$\$ where BGSS
is the between-group and WGSS the within-group sum of squares. The \\k\\
maximising CH is retained. Clusters are relabelled in increasing order
of their centres so that `Envirotype_<var>_00001` is always the lowest
range of the variable. For joint (multivariate) mining there is no
natural ordering, so centroids are ordered along their first principal
component and variables are standardised within the window before
clustering so that no feature dominates the Euclidean distance.

## See also

`W_matrix`, `env_kernel`,
[`T_matrix`](https://gcostaneto.github.io/EnvRtype/reference/T_matrix.md)

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
data("maizeWTH")
env.data <- maizeWTH[maizeWTH$daysFromStart < 100, ]

## 1. Generic time intervals
env_typing(env.data = env.data, env.id = 'env', var.id = 'T2M', by.interval = TRUE)
#> ---------------------------------------------------------------
#> env_typing -- mines environmental types from weather data
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assigning records to time-interval stages
#>                    env.variable Freq env    interval var
#> 1      T2M_(17.4,22]_Interval_0    0  NM  Interval_0 T2M
#> 2      T2M_(22,24.1]_Interval_0    3  NM  Interval_0 T2M
#> 3    T2M_(24.1,25.2]_Interval_0    4  NM  Interval_0 T2M
#> 4    T2M_(25.2,27.9]_Interval_0    2  NM  Interval_0 T2M
#> 5     T2M_(17.4,22]_Interval_10    0  NM Interval_10 T2M
#> 6     T2M_(22,24.1]_Interval_10    0  NM Interval_10 T2M
#> 7   T2M_(24.1,25.2]_Interval_10    8  NM Interval_10 T2M
#> 8   T2M_(25.2,27.9]_Interval_10    2  NM Interval_10 T2M
#> 9     T2M_(17.4,22]_Interval_20    0  NM Interval_20 T2M
#> 10    T2M_(22,24.1]_Interval_20    3  NM Interval_20 T2M
#> 11  T2M_(24.1,25.2]_Interval_20    5  NM Interval_20 T2M
#> 12  T2M_(25.2,27.9]_Interval_20    2  NM Interval_20 T2M
#> 13    T2M_(17.4,22]_Interval_30    0  NM Interval_30 T2M
#> 14    T2M_(22,24.1]_Interval_30    2  NM Interval_30 T2M
#> 15  T2M_(24.1,25.2]_Interval_30    7  NM Interval_30 T2M
#> 16  T2M_(25.2,27.9]_Interval_30    1  NM Interval_30 T2M
#> 17    T2M_(17.4,22]_Interval_40    0  NM Interval_40 T2M
#> 18    T2M_(22,24.1]_Interval_40    0  NM Interval_40 T2M
#> 19  T2M_(24.1,25.2]_Interval_40    4  NM Interval_40 T2M
#> 20  T2M_(25.2,27.9]_Interval_40    6  NM Interval_40 T2M
#> 21    T2M_(17.4,22]_Interval_50    0  NM Interval_50 T2M
#> 22    T2M_(22,24.1]_Interval_50    0  NM Interval_50 T2M
#> 23  T2M_(24.1,25.2]_Interval_50   10  NM Interval_50 T2M
#> 24  T2M_(25.2,27.9]_Interval_50    0  NM Interval_50 T2M
#> 25    T2M_(17.4,22]_Interval_60    0  NM Interval_60 T2M
#> 26    T2M_(22,24.1]_Interval_60    0  NM Interval_60 T2M
#> 27  T2M_(24.1,25.2]_Interval_60    7  NM Interval_60 T2M
#> 28  T2M_(25.2,27.9]_Interval_60    3  NM Interval_60 T2M
#> 29    T2M_(17.4,22]_Interval_70    0  NM Interval_70 T2M
#> 30    T2M_(22,24.1]_Interval_70    4  NM Interval_70 T2M
#> 31  T2M_(24.1,25.2]_Interval_70    6  NM Interval_70 T2M
#> 32  T2M_(25.2,27.9]_Interval_70    0  NM Interval_70 T2M
#> 33    T2M_(17.4,22]_Interval_80    0  NM Interval_80 T2M
#> 34    T2M_(22,24.1]_Interval_80    4  NM Interval_80 T2M
#> 35  T2M_(24.1,25.2]_Interval_80    6  NM Interval_80 T2M
#> 36  T2M_(25.2,27.9]_Interval_80    0  NM Interval_80 T2M
#> 37    T2M_(17.4,22]_Interval_90    0  NM Interval_90 T2M
#> 38    T2M_(22,24.1]_Interval_90    3  NM Interval_90 T2M
#> 39  T2M_(24.1,25.2]_Interval_90    6  NM Interval_90 T2M
#> 40  T2M_(25.2,27.9]_Interval_90    1  NM Interval_90 T2M
#> 41     T2M_(17.4,22]_Interval_0    0  PM  Interval_0 T2M
#> 42     T2M_(22,24.1]_Interval_0    7  PM  Interval_0 T2M
#> 43   T2M_(24.1,25.2]_Interval_0    2  PM  Interval_0 T2M
#> 44   T2M_(25.2,27.9]_Interval_0    0  PM  Interval_0 T2M
#> 45    T2M_(17.4,22]_Interval_10    5  PM Interval_10 T2M
#> 46    T2M_(22,24.1]_Interval_10    5  PM Interval_10 T2M
#> 47  T2M_(24.1,25.2]_Interval_10    0  PM Interval_10 T2M
#> 48  T2M_(25.2,27.9]_Interval_10    0  PM Interval_10 T2M
#> 49    T2M_(17.4,22]_Interval_20    4  PM Interval_20 T2M
#> 50    T2M_(22,24.1]_Interval_20    5  PM Interval_20 T2M
#> 51  T2M_(24.1,25.2]_Interval_20    1  PM Interval_20 T2M
#> 52  T2M_(25.2,27.9]_Interval_20    0  PM Interval_20 T2M
#> 53    T2M_(17.4,22]_Interval_30    1  PM Interval_30 T2M
#> 54    T2M_(22,24.1]_Interval_30    4  PM Interval_30 T2M
#> 55  T2M_(24.1,25.2]_Interval_30    5  PM Interval_30 T2M
#> 56  T2M_(25.2,27.9]_Interval_30    0  PM Interval_30 T2M
#> 57    T2M_(17.4,22]_Interval_40    3  PM Interval_40 T2M
#> 58    T2M_(22,24.1]_Interval_40    7  PM Interval_40 T2M
#> 59  T2M_(24.1,25.2]_Interval_40    0  PM Interval_40 T2M
#> 60  T2M_(25.2,27.9]_Interval_40    0  PM Interval_40 T2M
#> 61    T2M_(17.4,22]_Interval_50    1  PM Interval_50 T2M
#> 62    T2M_(22,24.1]_Interval_50    7  PM Interval_50 T2M
#> 63  T2M_(24.1,25.2]_Interval_50    2  PM Interval_50 T2M
#> 64  T2M_(25.2,27.9]_Interval_50    0  PM Interval_50 T2M
#> 65    T2M_(17.4,22]_Interval_60    9  PM Interval_60 T2M
#> 66    T2M_(22,24.1]_Interval_60    1  PM Interval_60 T2M
#> 67  T2M_(24.1,25.2]_Interval_60    0  PM Interval_60 T2M
#> 68  T2M_(25.2,27.9]_Interval_60    0  PM Interval_60 T2M
#> 69    T2M_(17.4,22]_Interval_70   10  PM Interval_70 T2M
#> 70    T2M_(22,24.1]_Interval_70    0  PM Interval_70 T2M
#> 71  T2M_(24.1,25.2]_Interval_70    0  PM Interval_70 T2M
#> 72  T2M_(25.2,27.9]_Interval_70    0  PM Interval_70 T2M
#> 73    T2M_(17.4,22]_Interval_80   10  PM Interval_80 T2M
#> 74    T2M_(22,24.1]_Interval_80    0  PM Interval_80 T2M
#> 75  T2M_(24.1,25.2]_Interval_80    0  PM Interval_80 T2M
#> 76  T2M_(25.2,27.9]_Interval_80    0  PM Interval_80 T2M
#> 77    T2M_(17.4,22]_Interval_90    9  PM Interval_90 T2M
#> 78    T2M_(22,24.1]_Interval_90    1  PM Interval_90 T2M
#> 79  T2M_(24.1,25.2]_Interval_90    0  PM Interval_90 T2M
#> 80  T2M_(25.2,27.9]_Interval_90    0  PM Interval_90 T2M
#> 81     T2M_(17.4,22]_Interval_0    0  IP  Interval_0 T2M
#> 82     T2M_(22,24.1]_Interval_0    3  IP  Interval_0 T2M
#> 83   T2M_(24.1,25.2]_Interval_0    3  IP  Interval_0 T2M
#> 84   T2M_(25.2,27.9]_Interval_0    3  IP  Interval_0 T2M
#> 85    T2M_(17.4,22]_Interval_10    0  IP Interval_10 T2M
#> 86    T2M_(22,24.1]_Interval_10    2  IP Interval_10 T2M
#> 87  T2M_(24.1,25.2]_Interval_10    6  IP Interval_10 T2M
#> 88  T2M_(25.2,27.9]_Interval_10    2  IP Interval_10 T2M
#> 89    T2M_(17.4,22]_Interval_20    0  IP Interval_20 T2M
#> 90    T2M_(22,24.1]_Interval_20    1  IP Interval_20 T2M
#> 91  T2M_(24.1,25.2]_Interval_20    7  IP Interval_20 T2M
#> 92  T2M_(25.2,27.9]_Interval_20    2  IP Interval_20 T2M
#> 93    T2M_(17.4,22]_Interval_30    1  IP Interval_30 T2M
#> 94    T2M_(22,24.1]_Interval_30    9  IP Interval_30 T2M
#> 95  T2M_(24.1,25.2]_Interval_30    0  IP Interval_30 T2M
#> 96  T2M_(25.2,27.9]_Interval_30    0  IP Interval_30 T2M
#> 97    T2M_(17.4,22]_Interval_40    0  IP Interval_40 T2M
#> 98    T2M_(22,24.1]_Interval_40    1  IP Interval_40 T2M
#> 99  T2M_(24.1,25.2]_Interval_40    5  IP Interval_40 T2M
#> 100 T2M_(25.2,27.9]_Interval_40    4  IP Interval_40 T2M
#> 101   T2M_(17.4,22]_Interval_50    0  IP Interval_50 T2M
#> 102   T2M_(22,24.1]_Interval_50    2  IP Interval_50 T2M
#> 103 T2M_(24.1,25.2]_Interval_50    7  IP Interval_50 T2M
#> 104 T2M_(25.2,27.9]_Interval_50    1  IP Interval_50 T2M
#> 105   T2M_(17.4,22]_Interval_60    0  IP Interval_60 T2M
#> 106   T2M_(22,24.1]_Interval_60    2  IP Interval_60 T2M
#> 107 T2M_(24.1,25.2]_Interval_60    6  IP Interval_60 T2M
#> 108 T2M_(25.2,27.9]_Interval_60    2  IP Interval_60 T2M
#> 109   T2M_(17.4,22]_Interval_70    0  IP Interval_70 T2M
#> 110   T2M_(22,24.1]_Interval_70    6  IP Interval_70 T2M
#> 111 T2M_(24.1,25.2]_Interval_70    3  IP Interval_70 T2M
#> 112 T2M_(25.2,27.9]_Interval_70    1  IP Interval_70 T2M
#> 113   T2M_(17.4,22]_Interval_80    2  IP Interval_80 T2M
#> 114   T2M_(22,24.1]_Interval_80    8  IP Interval_80 T2M
#> 115 T2M_(24.1,25.2]_Interval_80    0  IP Interval_80 T2M
#> 116 T2M_(25.2,27.9]_Interval_80    0  IP Interval_80 T2M
#> 117   T2M_(17.4,22]_Interval_90   10  IP Interval_90 T2M
#> 118   T2M_(22,24.1]_Interval_90    0  IP Interval_90 T2M
#> 119 T2M_(24.1,25.2]_Interval_90    0  IP Interval_90 T2M
#> 120 T2M_(25.2,27.9]_Interval_90    0  IP Interval_90 T2M
#> 121    T2M_(17.4,22]_Interval_0    0  SE  Interval_0 T2M
#> 122    T2M_(22,24.1]_Interval_0    6  SE  Interval_0 T2M
#> 123  T2M_(24.1,25.2]_Interval_0    3  SE  Interval_0 T2M
#> 124  T2M_(25.2,27.9]_Interval_0    0  SE  Interval_0 T2M
#> 125   T2M_(17.4,22]_Interval_10    3  SE Interval_10 T2M
#> 126   T2M_(22,24.1]_Interval_10    7  SE Interval_10 T2M
#> 127 T2M_(24.1,25.2]_Interval_10    0  SE Interval_10 T2M
#> 128 T2M_(25.2,27.9]_Interval_10    0  SE Interval_10 T2M
#> 129   T2M_(17.4,22]_Interval_20    3  SE Interval_20 T2M
#> 130   T2M_(22,24.1]_Interval_20    6  SE Interval_20 T2M
#> 131 T2M_(24.1,25.2]_Interval_20    1  SE Interval_20 T2M
#> 132 T2M_(25.2,27.9]_Interval_20    0  SE Interval_20 T2M
#> 133   T2M_(17.4,22]_Interval_30    6  SE Interval_30 T2M
#> 134   T2M_(22,24.1]_Interval_30    4  SE Interval_30 T2M
#> 135 T2M_(24.1,25.2]_Interval_30    0  SE Interval_30 T2M
#> 136 T2M_(25.2,27.9]_Interval_30    0  SE Interval_30 T2M
#> 137   T2M_(17.4,22]_Interval_40    0  SE Interval_40 T2M
#> 138   T2M_(22,24.1]_Interval_40    7  SE Interval_40 T2M
#> 139 T2M_(24.1,25.2]_Interval_40    3  SE Interval_40 T2M
#> 140 T2M_(25.2,27.9]_Interval_40    0  SE Interval_40 T2M
#> 141   T2M_(17.4,22]_Interval_50    9  SE Interval_50 T2M
#> 142   T2M_(22,24.1]_Interval_50    1  SE Interval_50 T2M
#> 143 T2M_(24.1,25.2]_Interval_50    0  SE Interval_50 T2M
#> 144 T2M_(25.2,27.9]_Interval_50    0  SE Interval_50 T2M
#> 145   T2M_(17.4,22]_Interval_60    8  SE Interval_60 T2M
#> 146   T2M_(22,24.1]_Interval_60    0  SE Interval_60 T2M
#> 147 T2M_(24.1,25.2]_Interval_60    0  SE Interval_60 T2M
#> 148 T2M_(25.2,27.9]_Interval_60    0  SE Interval_60 T2M
#> 149   T2M_(17.4,22]_Interval_70   10  SE Interval_70 T2M
#> 150   T2M_(22,24.1]_Interval_70    0  SE Interval_70 T2M
#> 151 T2M_(24.1,25.2]_Interval_70    0  SE Interval_70 T2M
#> 152 T2M_(25.2,27.9]_Interval_70    0  SE Interval_70 T2M
#> 153   T2M_(17.4,22]_Interval_80    6  SE Interval_80 T2M
#> 154   T2M_(22,24.1]_Interval_80    2  SE Interval_80 T2M
#> 155 T2M_(24.1,25.2]_Interval_80    0  SE Interval_80 T2M
#> 156 T2M_(25.2,27.9]_Interval_80    0  SE Interval_80 T2M
#> 157   T2M_(17.4,22]_Interval_90    9  SE Interval_90 T2M
#> 158   T2M_(22,24.1]_Interval_90    0  SE Interval_90 T2M
#> 159 T2M_(24.1,25.2]_Interval_90    0  SE Interval_90 T2M
#> 160 T2M_(25.2,27.9]_Interval_90    0  SE Interval_90 T2M
#> 161    T2M_(17.4,22]_Interval_0    0  SO  Interval_0 T2M
#> 162    T2M_(22,24.1]_Interval_0    0  SO  Interval_0 T2M
#> 163  T2M_(24.1,25.2]_Interval_0    0  SO  Interval_0 T2M
#> 164  T2M_(25.2,27.9]_Interval_0    9  SO  Interval_0 T2M
#> 165   T2M_(17.4,22]_Interval_10    0  SO Interval_10 T2M
#> 166   T2M_(22,24.1]_Interval_10    0  SO Interval_10 T2M
#> 167 T2M_(24.1,25.2]_Interval_10    0  SO Interval_10 T2M
#> 168 T2M_(25.2,27.9]_Interval_10    9  SO Interval_10 T2M
#> 169   T2M_(17.4,22]_Interval_20    0  SO Interval_20 T2M
#> 170   T2M_(22,24.1]_Interval_20    0  SO Interval_20 T2M
#> 171 T2M_(24.1,25.2]_Interval_20    0  SO Interval_20 T2M
#> 172 T2M_(25.2,27.9]_Interval_20    9  SO Interval_20 T2M
#> 173   T2M_(17.4,22]_Interval_30    0  SO Interval_30 T2M
#> 174   T2M_(22,24.1]_Interval_30    0  SO Interval_30 T2M
#> 175 T2M_(24.1,25.2]_Interval_30    0  SO Interval_30 T2M
#> 176 T2M_(25.2,27.9]_Interval_30    9  SO Interval_30 T2M
#> 177   T2M_(17.4,22]_Interval_40    0  SO Interval_40 T2M
#> 178   T2M_(22,24.1]_Interval_40    0  SO Interval_40 T2M
#> 179 T2M_(24.1,25.2]_Interval_40    2  SO Interval_40 T2M
#> 180 T2M_(25.2,27.9]_Interval_40    8  SO Interval_40 T2M
#> 181   T2M_(17.4,22]_Interval_50    0  SO Interval_50 T2M
#> 182   T2M_(22,24.1]_Interval_50    0  SO Interval_50 T2M
#> 183 T2M_(24.1,25.2]_Interval_50    0  SO Interval_50 T2M
#> 184 T2M_(25.2,27.9]_Interval_50   10  SO Interval_50 T2M
#> 185   T2M_(17.4,22]_Interval_60    0  SO Interval_60 T2M
#> 186   T2M_(22,24.1]_Interval_60    1  SO Interval_60 T2M
#> 187 T2M_(24.1,25.2]_Interval_60    3  SO Interval_60 T2M
#> 188 T2M_(25.2,27.9]_Interval_60    6  SO Interval_60 T2M
#> 189   T2M_(17.4,22]_Interval_70    0  SO Interval_70 T2M
#> 190   T2M_(22,24.1]_Interval_70    0  SO Interval_70 T2M
#> 191 T2M_(24.1,25.2]_Interval_70    2  SO Interval_70 T2M
#> 192 T2M_(25.2,27.9]_Interval_70    8  SO Interval_70 T2M
#> 193   T2M_(17.4,22]_Interval_80    0  SO Interval_80 T2M
#> 194   T2M_(22,24.1]_Interval_80    0  SO Interval_80 T2M
#> 195 T2M_(24.1,25.2]_Interval_80    0  SO Interval_80 T2M
#> 196 T2M_(25.2,27.9]_Interval_80   10  SO Interval_80 T2M
#> 197   T2M_(17.4,22]_Interval_90    0  SO Interval_90 T2M
#> 198   T2M_(22,24.1]_Interval_90    0  SO Interval_90 T2M
#> 199 T2M_(24.1,25.2]_Interval_90    0  SO Interval_90 T2M
#> 200 T2M_(25.2,27.9]_Interval_90    8  SO Interval_90 T2M

## 2. Specific time intervals with names
env_typing(env.data = env.data, env.id = 'env', var.id = 'T2M',
           by.interval = TRUE,
           time.window  = c(0, 15, 35, 65, 90),
           names.window = c('1-initial growing', '2-leaf expansion I',
                            '3-leaf expansion II', '4-flowering', '5-grain filling'))
#> ---------------------------------------------------------------
#> env_typing -- mines environmental types from weather data
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assigning records to time-interval stages
#>                            env.variable Freq env            interval var
#> 1       T2M_(17.4,22]_1-initial growing    0  NM   1-initial growing T2M
#> 2       T2M_(22,24.1]_1-initial growing    3  NM   1-initial growing T2M
#> 3     T2M_(24.1,25.2]_1-initial growing    8  NM   1-initial growing T2M
#> 4     T2M_(25.2,27.9]_1-initial growing    3  NM   1-initial growing T2M
#> 5      T2M_(17.4,22]_2-leaf expansion I    0  NM  2-leaf expansion I T2M
#> 6      T2M_(22,24.1]_2-leaf expansion I    5  NM  2-leaf expansion I T2M
#> 7    T2M_(24.1,25.2]_2-leaf expansion I   12  NM  2-leaf expansion I T2M
#> 8    T2M_(25.2,27.9]_2-leaf expansion I    3  NM  2-leaf expansion I T2M
#> 9     T2M_(17.4,22]_3-leaf expansion II    0  NM 3-leaf expansion II T2M
#> 10    T2M_(22,24.1]_3-leaf expansion II    0  NM 3-leaf expansion II T2M
#> 11  T2M_(24.1,25.2]_3-leaf expansion II   22  NM 3-leaf expansion II T2M
#> 12  T2M_(25.2,27.9]_3-leaf expansion II    8  NM 3-leaf expansion II T2M
#> 13            T2M_(17.4,22]_4-flowering    0  NM         4-flowering T2M
#> 14            T2M_(22,24.1]_4-flowering    8  NM         4-flowering T2M
#> 15          T2M_(24.1,25.2]_4-flowering   15  NM         4-flowering T2M
#> 16          T2M_(25.2,27.9]_4-flowering    2  NM         4-flowering T2M
#> 17        T2M_(17.4,22]_5-grain filling    0  NM     5-grain filling T2M
#> 18        T2M_(22,24.1]_5-grain filling    3  NM     5-grain filling T2M
#> 19      T2M_(24.1,25.2]_5-grain filling    6  NM     5-grain filling T2M
#> 20      T2M_(25.2,27.9]_5-grain filling    1  NM     5-grain filling T2M
#> 21      T2M_(17.4,22]_1-initial growing    2  PM   1-initial growing T2M
#> 22      T2M_(22,24.1]_1-initial growing   10  PM   1-initial growing T2M
#> 23    T2M_(24.1,25.2]_1-initial growing    2  PM   1-initial growing T2M
#> 24    T2M_(25.2,27.9]_1-initial growing    0  PM   1-initial growing T2M
#> 25     T2M_(17.4,22]_2-leaf expansion I    7  PM  2-leaf expansion I T2M
#> 26     T2M_(22,24.1]_2-leaf expansion I    7  PM  2-leaf expansion I T2M
#> 27   T2M_(24.1,25.2]_2-leaf expansion I    6  PM  2-leaf expansion I T2M
#> 28   T2M_(25.2,27.9]_2-leaf expansion I    0  PM  2-leaf expansion I T2M
#> 29    T2M_(17.4,22]_3-leaf expansion II   10  PM 3-leaf expansion II T2M
#> 30    T2M_(22,24.1]_3-leaf expansion II   18  PM 3-leaf expansion II T2M
#> 31  T2M_(24.1,25.2]_3-leaf expansion II    2  PM 3-leaf expansion II T2M
#> 32  T2M_(25.2,27.9]_3-leaf expansion II    0  PM 3-leaf expansion II T2M
#> 33            T2M_(17.4,22]_4-flowering   24  PM         4-flowering T2M
#> 34            T2M_(22,24.1]_4-flowering    1  PM         4-flowering T2M
#> 35          T2M_(24.1,25.2]_4-flowering    0  PM         4-flowering T2M
#> 36          T2M_(25.2,27.9]_4-flowering    0  PM         4-flowering T2M
#> 37        T2M_(17.4,22]_5-grain filling    9  PM     5-grain filling T2M
#> 38        T2M_(22,24.1]_5-grain filling    1  PM     5-grain filling T2M
#> 39      T2M_(24.1,25.2]_5-grain filling    0  PM     5-grain filling T2M
#> 40      T2M_(25.2,27.9]_5-grain filling    0  PM     5-grain filling T2M
#> 41      T2M_(17.4,22]_1-initial growing    0  IP   1-initial growing T2M
#> 42      T2M_(22,24.1]_1-initial growing    5  IP   1-initial growing T2M
#> 43    T2M_(24.1,25.2]_1-initial growing    5  IP   1-initial growing T2M
#> 44    T2M_(25.2,27.9]_1-initial growing    4  IP   1-initial growing T2M
#> 45     T2M_(17.4,22]_2-leaf expansion I    0  IP  2-leaf expansion I T2M
#> 46     T2M_(22,24.1]_2-leaf expansion I    6  IP  2-leaf expansion I T2M
#> 47   T2M_(24.1,25.2]_2-leaf expansion I   11  IP  2-leaf expansion I T2M
#> 48   T2M_(25.2,27.9]_2-leaf expansion I    3  IP  2-leaf expansion I T2M
#> 49    T2M_(17.4,22]_3-leaf expansion II    1  IP 3-leaf expansion II T2M
#> 50    T2M_(22,24.1]_3-leaf expansion II    9  IP 3-leaf expansion II T2M
#> 51  T2M_(24.1,25.2]_3-leaf expansion II   15  IP 3-leaf expansion II T2M
#> 52  T2M_(25.2,27.9]_3-leaf expansion II    5  IP 3-leaf expansion II T2M
#> 53            T2M_(17.4,22]_4-flowering    2  IP         4-flowering T2M
#> 54            T2M_(22,24.1]_4-flowering   14  IP         4-flowering T2M
#> 55          T2M_(24.1,25.2]_4-flowering    6  IP         4-flowering T2M
#> 56          T2M_(25.2,27.9]_4-flowering    3  IP         4-flowering T2M
#> 57        T2M_(17.4,22]_5-grain filling   10  IP     5-grain filling T2M
#> 58        T2M_(22,24.1]_5-grain filling    0  IP     5-grain filling T2M
#> 59      T2M_(24.1,25.2]_5-grain filling    0  IP     5-grain filling T2M
#> 60      T2M_(25.2,27.9]_5-grain filling    0  IP     5-grain filling T2M
#> 61      T2M_(17.4,22]_1-initial growing    0  SE   1-initial growing T2M
#> 62      T2M_(22,24.1]_1-initial growing   11  SE   1-initial growing T2M
#> 63    T2M_(24.1,25.2]_1-initial growing    3  SE   1-initial growing T2M
#> 64    T2M_(25.2,27.9]_1-initial growing    0  SE   1-initial growing T2M
#> 65     T2M_(17.4,22]_2-leaf expansion I    9  SE  2-leaf expansion I T2M
#> 66     T2M_(22,24.1]_2-leaf expansion I   10  SE  2-leaf expansion I T2M
#> 67   T2M_(24.1,25.2]_2-leaf expansion I    1  SE  2-leaf expansion I T2M
#> 68   T2M_(25.2,27.9]_2-leaf expansion I    0  SE  2-leaf expansion I T2M
#> 69    T2M_(17.4,22]_3-leaf expansion II   17  SE 3-leaf expansion II T2M
#> 70    T2M_(22,24.1]_3-leaf expansion II   10  SE 3-leaf expansion II T2M
#> 71  T2M_(24.1,25.2]_3-leaf expansion II    3  SE 3-leaf expansion II T2M
#> 72  T2M_(25.2,27.9]_3-leaf expansion II    0  SE 3-leaf expansion II T2M
#> 73            T2M_(17.4,22]_4-flowering   19  SE         4-flowering T2M
#> 74            T2M_(22,24.1]_4-flowering    2  SE         4-flowering T2M
#> 75          T2M_(24.1,25.2]_4-flowering    0  SE         4-flowering T2M
#> 76          T2M_(25.2,27.9]_4-flowering    0  SE         4-flowering T2M
#> 77        T2M_(17.4,22]_5-grain filling    9  SE     5-grain filling T2M
#> 78        T2M_(22,24.1]_5-grain filling    0  SE     5-grain filling T2M
#> 79      T2M_(24.1,25.2]_5-grain filling    0  SE     5-grain filling T2M
#> 80      T2M_(25.2,27.9]_5-grain filling    0  SE     5-grain filling T2M
#> 81      T2M_(17.4,22]_1-initial growing    0  SO   1-initial growing T2M
#> 82      T2M_(22,24.1]_1-initial growing    0  SO   1-initial growing T2M
#> 83    T2M_(24.1,25.2]_1-initial growing    0  SO   1-initial growing T2M
#> 84    T2M_(25.2,27.9]_1-initial growing   13  SO   1-initial growing T2M
#> 85     T2M_(17.4,22]_2-leaf expansion I    0  SO  2-leaf expansion I T2M
#> 86     T2M_(22,24.1]_2-leaf expansion I    0  SO  2-leaf expansion I T2M
#> 87   T2M_(24.1,25.2]_2-leaf expansion I    0  SO  2-leaf expansion I T2M
#> 88   T2M_(25.2,27.9]_2-leaf expansion I   18  SO  2-leaf expansion I T2M
#> 89    T2M_(17.4,22]_3-leaf expansion II    0  SO 3-leaf expansion II T2M
#> 90    T2M_(22,24.1]_3-leaf expansion II    0  SO 3-leaf expansion II T2M
#> 91  T2M_(24.1,25.2]_3-leaf expansion II    3  SO 3-leaf expansion II T2M
#> 92  T2M_(25.2,27.9]_3-leaf expansion II   27  SO 3-leaf expansion II T2M
#> 93            T2M_(17.4,22]_4-flowering    0  SO         4-flowering T2M
#> 94            T2M_(22,24.1]_4-flowering    1  SO         4-flowering T2M
#> 95          T2M_(24.1,25.2]_4-flowering    4  SO         4-flowering T2M
#> 96          T2M_(25.2,27.9]_4-flowering   20  SO         4-flowering T2M
#> 97        T2M_(17.4,22]_5-grain filling    0  SO     5-grain filling T2M
#> 98        T2M_(22,24.1]_5-grain filling    0  SO     5-grain filling T2M
#> 99      T2M_(24.1,25.2]_5-grain filling    0  SO     5-grain filling T2M
#> 100     T2M_(25.2,27.9]_5-grain filling    8  SO     5-grain filling T2M

## 3. Cardinal (ecophysiological) thresholds
env_typing(env.data = env.data, env.id = 'env',
           var.id    = c('T2M', 'PRECTOT'),
           cardinals = list(T2M = c(0, 9, 22, 32, 45), PRECTOT = c(0, 5, 10, 25, 100)))
#> ---------------------------------------------------------------
#> env_typing -- mines environmental types from weather data
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assigning records to time-interval stages
#>        env.variable Freq env       interval     var
#> 1         T2M_(0,9]    0  NM by environment     T2M
#> 2        T2M_(9,22]    0  NM by environment     T2M
#> 3       T2M_(22,32]   99  NM by environment     T2M
#> 4       T2M_(32,45]    0  NM by environment     T2M
#> 5         T2M_(0,9]    0  PM by environment     T2M
#> 6        T2M_(9,22]   53  PM by environment     T2M
#> 7       T2M_(22,32]   46  PM by environment     T2M
#> 8       T2M_(32,45]    0  PM by environment     T2M
#> 9         T2M_(0,9]    0  IP by environment     T2M
#> 10       T2M_(9,22]   13  IP by environment     T2M
#> 11      T2M_(22,32]   86  IP by environment     T2M
#> 12      T2M_(32,45]    0  IP by environment     T2M
#> 13        T2M_(0,9]    0  SE by environment     T2M
#> 14       T2M_(9,22]   60  SE by environment     T2M
#> 15      T2M_(22,32]   39  SE by environment     T2M
#> 16      T2M_(32,45]    0  SE by environment     T2M
#> 17        T2M_(0,9]    0  SO by environment     T2M
#> 18       T2M_(9,22]    0  SO by environment     T2M
#> 19      T2M_(22,32]   99  SO by environment     T2M
#> 20      T2M_(32,45]    0  SO by environment     T2M
#> 21    PRECTOT_(0,5]   70  NM by environment PRECTOT
#> 22   PRECTOT_(5,10]   14  NM by environment PRECTOT
#> 23  PRECTOT_(10,25]    8  NM by environment PRECTOT
#> 24 PRECTOT_(25,100]    1  NM by environment PRECTOT
#> 25    PRECTOT_(0,5]   62  PM by environment PRECTOT
#> 26   PRECTOT_(5,10]   14  PM by environment PRECTOT
#> 27  PRECTOT_(10,25]    6  PM by environment PRECTOT
#> 28 PRECTOT_(25,100]    0  PM by environment PRECTOT
#> 29    PRECTOT_(0,5]   68  IP by environment PRECTOT
#> 30   PRECTOT_(5,10]   16  IP by environment PRECTOT
#> 31  PRECTOT_(10,25]   13  IP by environment PRECTOT
#> 32 PRECTOT_(25,100]    0  IP by environment PRECTOT
#> 33    PRECTOT_(0,5]   53  SE by environment PRECTOT
#> 34   PRECTOT_(5,10]   12  SE by environment PRECTOT
#> 35  PRECTOT_(10,25]    8  SE by environment PRECTOT
#> 36 PRECTOT_(25,100]    1  SE by environment PRECTOT
#> 37    PRECTOT_(0,5]   80  SO by environment PRECTOT
#> 38   PRECTOT_(5,10]    8  SO by environment PRECTOT
#> 39  PRECTOT_(10,25]    2  SO by environment PRECTOT
#> 40 PRECTOT_(25,100]    0  SO by environment PRECTOT

## 4. Data-driven univariate mining (K-means + Calinski-Harabasz)
out <- env_typing(env.data = env.data, env.id = 'env',
                  var.id = c('T2M', 'PRECTOT'),
                  envirotype_mining = TRUE, k.range = 2:8)
#> ---------------------------------------------------------------
#> env_typing -- mines environmental types from weather data
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assigning records to time-interval stages
#> Mining envirotypes for: T2M
#> Mining envirotypes for: PRECTOT
out$typologies
#>                env.variable Freq env       interval     var
#> 1      Envirotype_T2M_00001    0  NM by environment     T2M
#> 2      Envirotype_T2M_00002    0  NM by environment     T2M
#> 3      Envirotype_T2M_00003    0  NM by environment     T2M
#> 4      Envirotype_T2M_00004    8  NM by environment     T2M
#> 5      Envirotype_T2M_00005   37  NM by environment     T2M
#> 6      Envirotype_T2M_00006   49  NM by environment     T2M
#> 7      Envirotype_T2M_00007    5  NM by environment     T2M
#> 8      Envirotype_T2M_00008    0  NM by environment     T2M
#> 9      Envirotype_T2M_00001    0  PM by environment     T2M
#> 10     Envirotype_T2M_00002   17  PM by environment     T2M
#> 11     Envirotype_T2M_00003   37  PM by environment     T2M
#> 12     Envirotype_T2M_00004   23  PM by environment     T2M
#> 13     Envirotype_T2M_00005   19  PM by environment     T2M
#> 14     Envirotype_T2M_00006    3  PM by environment     T2M
#> 15     Envirotype_T2M_00007    0  PM by environment     T2M
#> 16     Envirotype_T2M_00008    0  PM by environment     T2M
#> 17     Envirotype_T2M_00001    0  IP by environment     T2M
#> 18     Envirotype_T2M_00002    2  IP by environment     T2M
#> 19     Envirotype_T2M_00003   12  IP by environment     T2M
#> 20     Envirotype_T2M_00004   22  IP by environment     T2M
#> 21     Envirotype_T2M_00005   20  IP by environment     T2M
#> 22     Envirotype_T2M_00006   34  IP by environment     T2M
#> 23     Envirotype_T2M_00007    8  IP by environment     T2M
#> 24     Envirotype_T2M_00008    1  IP by environment     T2M
#> 25     Envirotype_T2M_00001   14  SE by environment     T2M
#> 26     Envirotype_T2M_00002   25  SE by environment     T2M
#> 27     Envirotype_T2M_00003   22  SE by environment     T2M
#> 28     Envirotype_T2M_00004   24  SE by environment     T2M
#> 29     Envirotype_T2M_00005   10  SE by environment     T2M
#> 30     Envirotype_T2M_00006    4  SE by environment     T2M
#> 31     Envirotype_T2M_00007    0  SE by environment     T2M
#> 32     Envirotype_T2M_00008    0  SE by environment     T2M
#> 33     Envirotype_T2M_00001    0  SO by environment     T2M
#> 34     Envirotype_T2M_00002    0  SO by environment     T2M
#> 35     Envirotype_T2M_00003    0  SO by environment     T2M
#> 36     Envirotype_T2M_00004    1  SO by environment     T2M
#> 37     Envirotype_T2M_00005    1  SO by environment     T2M
#> 38     Envirotype_T2M_00006   12  SO by environment     T2M
#> 39     Envirotype_T2M_00007   44  SO by environment     T2M
#> 40     Envirotype_T2M_00008   41  SO by environment     T2M
#> 41 Envirotype_PRECTOT_00001   39  NM by environment PRECTOT
#> 42 Envirotype_PRECTOT_00002   24  NM by environment PRECTOT
#> 43 Envirotype_PRECTOT_00003   15  NM by environment PRECTOT
#> 44 Envirotype_PRECTOT_00004   12  NM by environment PRECTOT
#> 45 Envirotype_PRECTOT_00005    4  NM by environment PRECTOT
#> 46 Envirotype_PRECTOT_00006    4  NM by environment PRECTOT
#> 47 Envirotype_PRECTOT_00007    0  NM by environment PRECTOT
#> 48 Envirotype_PRECTOT_00008    1  NM by environment PRECTOT
#> 49 Envirotype_PRECTOT_00001   51  PM by environment PRECTOT
#> 50 Envirotype_PRECTOT_00002   20  PM by environment PRECTOT
#> 51 Envirotype_PRECTOT_00003   11  PM by environment PRECTOT
#> 52 Envirotype_PRECTOT_00004    8  PM by environment PRECTOT
#> 53 Envirotype_PRECTOT_00005    8  PM by environment PRECTOT
#> 54 Envirotype_PRECTOT_00006    1  PM by environment PRECTOT
#> 55 Envirotype_PRECTOT_00007    0  PM by environment PRECTOT
#> 56 Envirotype_PRECTOT_00008    0  PM by environment PRECTOT
#> 57 Envirotype_PRECTOT_00001   37  IP by environment PRECTOT
#> 58 Envirotype_PRECTOT_00002   24  IP by environment PRECTOT
#> 59 Envirotype_PRECTOT_00003   12  IP by environment PRECTOT
#> 60 Envirotype_PRECTOT_00004   11  IP by environment PRECTOT
#> 61 Envirotype_PRECTOT_00005    6  IP by environment PRECTOT
#> 62 Envirotype_PRECTOT_00006    5  IP by environment PRECTOT
#> 63 Envirotype_PRECTOT_00007    4  IP by environment PRECTOT
#> 64 Envirotype_PRECTOT_00008    0  IP by environment PRECTOT
#> 65 Envirotype_PRECTOT_00001   64  SE by environment PRECTOT
#> 66 Envirotype_PRECTOT_00002    9  SE by environment PRECTOT
#> 67 Envirotype_PRECTOT_00003    6  SE by environment PRECTOT
#> 68 Envirotype_PRECTOT_00004   10  SE by environment PRECTOT
#> 69 Envirotype_PRECTOT_00005    4  SE by environment PRECTOT
#> 70 Envirotype_PRECTOT_00006    2  SE by environment PRECTOT
#> 71 Envirotype_PRECTOT_00007    3  SE by environment PRECTOT
#> 72 Envirotype_PRECTOT_00008    1  SE by environment PRECTOT
#> 73 Envirotype_PRECTOT_00001   71  SO by environment PRECTOT
#> 74 Envirotype_PRECTOT_00002   13  SO by environment PRECTOT
#> 75 Envirotype_PRECTOT_00003    6  SO by environment PRECTOT
#> 76 Envirotype_PRECTOT_00004    6  SO by environment PRECTOT
#> 77 Envirotype_PRECTOT_00005    2  SO by environment PRECTOT
#> 78 Envirotype_PRECTOT_00006    0  SO by environment PRECTOT
#> 79 Envirotype_PRECTOT_00007    1  SO by environment PRECTOT
#> 80 Envirotype_PRECTOT_00008    0  SO by environment PRECTOT
out$envirotype_description
#>    variable               envirotype cluster   n   min   max       mean
#> 1       T2M     Envirotype_T2M_00001       1  14 15.09 18.46 17.4542857
#> 2       T2M     Envirotype_T2M_00002       2  44 18.65 20.52 19.6634091
#> 3       T2M     Envirotype_T2M_00003       3  71 20.60 22.09 21.4129577
#> 4       T2M     Envirotype_T2M_00004       4  78 22.16 23.41 22.8153846
#> 5       T2M     Envirotype_T2M_00005       5  87 23.44 24.49 24.0488506
#> 6       T2M     Envirotype_T2M_00006       6 102 24.52 25.49 24.9660784
#> 7       T2M     Envirotype_T2M_00007       7  57 25.50 26.59 26.0212281
#> 8       T2M     Envirotype_T2M_00008       8  42 26.73 28.56 27.4154762
#> 9   PRECTOT Envirotype_PRECTOT_00001       1 262  0.00  1.35  0.3102672
#> 10  PRECTOT Envirotype_PRECTOT_00002       2  90  1.39  3.40  2.4115556
#> 11  PRECTOT Envirotype_PRECTOT_00003       3  50  3.46  5.50  4.4398000
#> 12  PRECTOT Envirotype_PRECTOT_00004       4  47  5.76  8.65  7.0387234
#> 13  PRECTOT Envirotype_PRECTOT_00005       5  24  9.07 12.89 10.8316667
#> 14  PRECTOT Envirotype_PRECTOT_00006       6  12 13.28 18.43 15.4133333
#> 15  PRECTOT Envirotype_PRECTOT_00007       7   8 19.52 24.92 23.0287500
#> 16  PRECTOT Envirotype_PRECTOT_00008       8   2 41.67 42.90 42.2850000
#>           sd median     center  lower  upper
#> 1  0.9307294 17.605 17.4542857   -Inf 18.555
#> 2  0.5449683 19.725 19.6634091 18.555 20.560
#> 3  0.3983516 21.490 21.4129577 20.560 22.125
#> 4  0.3728643 22.865 22.8153846 22.125 23.425
#> 5  0.2857525 24.090 24.0488506 23.425 24.505
#> 6  0.2470425 24.945 24.9660784 24.505 25.495
#> 7  0.3540534 25.970 26.0212281 25.495 26.660
#> 8  0.4456358 27.355 27.4154762 26.660    Inf
#> 9  0.3746195  0.160  0.3102672   -Inf  1.370
#> 10 0.6606056  2.460  2.4115556  1.370  3.430
#> 11 0.5705976  4.460  4.4398000  3.430  5.630
#> 12 0.8920412  6.960  7.0387234  5.630  8.860
#> 13 1.1609579 10.810 10.8316667  8.860 13.085
#> 14 1.3932065 15.250 15.4133333 13.085 18.975
#> 15 1.8954414 23.595 23.0287500 18.975 33.295
#> 16 0.8697413 42.285 42.2850000 33.295    Inf
out$mining_summary
#>   variable k.optimal CH.index   k.evaluated
#> 1      T2M         8 2580.297 2,3,4,5,6,7,8
#> 2  PRECTOT         8 4182.403 2,3,4,5,6,7,8
#>                                                       CH.profile
#> 1 1004.275,1063.028,1480.092,1731.761,1862.231,2246.823,2580.297
#> 2   790.549,1038.29,1324.786,1427.687,1413.491,3407.667,4182.403

## 5. Joint (multi-variable) envirotypes per development stage
outj <- env_typing(env.data = env.data, env.id = 'env',
                   var.id = c('T2M', 'PRECTOT'),
                   by.interval = TRUE, time.window = c(0, 30, 60, 90),
                   envirotype_mining = TRUE, joint = TRUE)
#> ---------------------------------------------------------------
#> env_typing -- mines environmental types from weather data
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assigning records to time-interval stages
#> Mining joint envirotypes for window: Interval_0
#> Mining joint envirotypes for window: Interval_30
#> Mining joint envirotypes for window: Interval_60
#> Mining joint envirotypes for window: Interval_90
head(outj$envirotype_description)
#>       window                  envirotype cluster variable  n   min   max
#> 1 Interval_0 Envirotype_Interval_0_00001       1      T2M 20 26.26 28.56
#> 2 Interval_0 Envirotype_Interval_0_00001       1  PRECTOT 20  0.01  1.06
#> 3 Interval_0 Envirotype_Interval_0_00002       2      T2M 25 24.49 25.97
#> 4 Interval_0 Envirotype_Interval_0_00002       2  PRECTOT 25  0.00  3.50
#> 5 Interval_0 Envirotype_Interval_0_00003       3      T2M 21 22.86 24.33
#> 6 Interval_0 Envirotype_Interval_0_00003       3  PRECTOT 21  0.20  4.15
#>        mean        sd median    center
#> 1 26.924000 0.6488970 26.655 26.924000
#> 2  0.463500 0.3212521  0.375  0.463500
#> 3 25.134000 0.4637887 25.050 25.134000
#> 4  1.273200 1.1913988  0.610  1.273200
#> 5 23.601429 0.4693110 23.540 23.601429
#> 6  1.472381 1.1584986  1.300  1.472381

## 6. Relative frequencies, wide format (ready for env_kernel)
ET <- env_typing(env.data = env.data, env.id = 'env', var.id = 'T2M',
                 format = 'wide', ratio = TRUE)
#> ---------------------------------------------------------------
#> env_typing -- mines environmental types from weather data
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assigning records to time-interval stages

## 7. Envirotypes for feature combinations (singles + pairs by default)
outc <- env_typing(env.data = env.data, env.id = 'env',
                   var.id = c('T2M', 'PRECTOT', 'VPD'),
                   envirotype_mining = TRUE, combine.features = TRUE)
#> ---------------------------------------------------------------
#> env_typing -- mines environmental types from weather data
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assigning records to time-interval stages
#> Mining combination FC001 (T2M)
#> Mining combination FC002 (PRECTOT)
#> Mining combination FC003 (VPD)
#> Mining combination FC004 (T2M+PRECTOT)
#> Mining combination FC005 (T2M+VPD)
#> Mining combination FC006 (PRECTOT+VPD)
outc$feature_combinations        # generic id -> which variables
#>   combo_id n_features     features      intervals
#> 1    FC001          1          T2M by environment
#> 2    FC002          1      PRECTOT by environment
#> 3    FC003          1          VPD by environment
#> 4    FC004          2 T2M, PRECTOT by environment
#> 5    FC005          2     T2M, VPD by environment
#> 6    FC006          2 PRECTOT, VPD by environment
head(outc$envirotype_description)
#>   combo_id features             envirotype cluster variable  n   min   max
#> 1    FC001      T2M Envirotype_FC001_00001       1      T2M  8 15.09 17.64
#> 2    FC001      T2M Envirotype_FC001_00002       2      T2M 23 17.96 19.48
#> 3    FC001      T2M Envirotype_FC001_00003       3      T2M 33 19.55 20.74
#> 4    FC001      T2M Envirotype_FC001_00004       4      T2M 55 20.80 21.84
#> 5    FC001      T2M Envirotype_FC001_00005       5      T2M 40 21.86 22.69
#> 6    FC001      T2M Envirotype_FC001_00006       6      T2M 57 22.73 23.62
#>       mean        sd median   center
#> 1 16.88625 0.8502426  17.12 16.88625
#> 2 18.84696 0.4471570  18.91 18.84696
#> 3 20.15091 0.3567944  20.11 20.15091
#> 4 21.39436 0.2936789  21.46 21.39436
#> 5 22.29625 0.2336301  22.31 22.29625
#> 6 23.14105 0.2513868  23.13 23.14105

## include triples too (combine.max.order = 3), or NULL for all orders
outc3 <- env_typing(env.data = env.data, env.id = 'env',
                    var.id = c('T2M', 'PRECTOT', 'VPD'),
                    envirotype_mining = TRUE, combine.features = TRUE,
                    combine.max.order = 3)
#> ---------------------------------------------------------------
#> env_typing -- mines environmental types from weather data
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assigning records to time-interval stages
#> Mining combination FC001 (T2M)
#> Mining combination FC002 (PRECTOT)
#> Mining combination FC003 (VPD)
#> Mining combination FC004 (T2M+PRECTOT)
#> Mining combination FC005 (T2M+VPD)
#> Mining combination FC006 (PRECTOT+VPD)
#> Mining combination FC007 (T2M+PRECTOT+VPD)
# }
```
