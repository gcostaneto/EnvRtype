# Build Interval-Resolved Environmental Indices

Converts daily environmental data into the environment x index matrix
used by Della Coletta et al. (2023): each weather factor is summarised
within fixed-width time intervals from planting, producing one column
per *factor x interval* combination.

## Usage

``` r
env_indices(
  env.data,
  env.id = "env",
  day.id = "daysFromStart",
  var.id = NULL,
  interval = 3L,
  end.day = 151L,
  statistic = c("mean", "sum", "min", "max", "median"),
  statistic.by = NULL,
  id.names = NULL,
  sep = "_",
  require.complete = TRUE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame of daily weather, e.g. a
  [`get_weather`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md)
  or
  [`processWTH`](https://gcostaneto.github.io/EnvRtype/reference/processWTH.md)
  output.

- env.id:

  character. Column identifying the environment. Default `"env"`.

- day.id:

  character. Column giving days from planting. Default
  `"daysFromStart"`.

- var.id:

  character vector. Which weather factors to use. `NULL` (default) uses
  every numeric column that is not an id column.

- interval:

  integer. Width of each time window in days. Default 3, as in the
  paper.

- end.day:

  integer. Last day retained. Default 151, as in the paper. Days beyond
  this are discarded so that every environment contributes the same
  intervals.

- statistic:

  character. Within-interval summary: `"mean"` (default), `"sum"`,
  `"min"`, `"max"` or `"median"`. Precipitation is arguably better
  summed than averaged – see `statistic.by`.

- statistic.by:

  named character vector. Per-factor override of `statistic`, e.g.
  `c(PRECTOT = "sum")`. Names are factor names, values are statistics.

- id.names:

  character vector. Additional columns to treat as identifiers (never
  summarised). Default `NULL`.

- sep:

  character. Separator between factor and interval in the output column
  names. Default `"_"`.

- require.complete:

  boolean. If `TRUE` (default) an error is raised when any environment x
  interval x factor cell is missing, because an incomplete matrix
  silently biases the PCA. Set `FALSE` to allow `NA`s through.

- verbose:

  boolean. Print a summary. Default `TRUE`.

## Value

A numeric matrix with environments in rows and \\n\_{factor} \times
n\_{interval}\\ indices in columns. Attributes:

- `"factors"`:

  the factor names used

- `"intervals"`:

  the interval labels

- `"interval.mid"`:

  named numeric vector of interval mid-point days

- `"interval.range"`:

  data.frame of interval start/end days

- `"n.missing"`:

  number of missing cells

## Details

**Why a fixed end day.** Environments differ in season length. Without a
common cut-off, environments contribute different numbers of intervals
and the matrix is ragged. The paper fixed the last interval at 151 days
after planting; the same device is used here via `end.day`.

**Interval labelling.** With `interval = 3` and `end.day = 151`, days
1-3 form interval 1, days 4-6 interval 2, and so on; the final partial
window (day 151) forms interval 51. This reproduces the 51 intervals of
the paper, and with 17 factors gives the 867 indices reported.

**Mean vs sum.** The paper does not state the within-interval statistic.
A mean is used by default because it is defined for every factor, but
for accumulating quantities such as precipitation a sum is more
interpretable; `statistic.by` allows a per-factor choice.

## References

Della Coletta R., Liese S.E., Fernandes S.B., Mikel M.A., Bohn M.O.,
Lipka A.E., Hirsch C.N. (2023). Linking genetic and environmental
factors through marker effect networks to understand trait plasticity.
*GENETICS* 224(4), iyad103.

## See also

[`env_pca`](https://gcostaneto.github.io/EnvRtype/reference/env_pca.md),
[`env_cor_heatmap`](https://gcostaneto.github.io/EnvRtype/reference/env_cor_heatmap.md)

## Author

Implementation following Della Coletta et al. (2023)

## Examples

``` r
# \donttest{
data("maizeWTH")

## 3-day intervals across the shared season. The bundled sample does not
## reach the paper's 151 days in every environment, so end.day is trimmed to
## a fully observed window (the default 151 would leave gaps and error out).
W <- env_indices(maizeWTH, env.id = "env", day.id = "daysFromStart",
                 end.day = 130)
#> ---------------------------------------------------------------
#> env_indices -- builds interval-resolved environmental indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ------------------------------------------------------------
#> env_indices(): 5 environments x 21 factors x 44 intervals = 924 indices
#>   interval width : 3 day(s), days 1..130
#>   discarded      : 106 record(s) outside 1..130
#> ------------------------------------------------------------
dim(W)
#> [1]   5 924
attr(W, "intervals")[1:5]
#> [1] "i001" "i002" "i003" "i004" "i005"

## Sum precipitation within each window, average everything else
W2 <- env_indices(maizeWTH, statistic = "mean",
                  statistic.by = c(PRECTOT = "sum"), end.day = 130)
#> ---------------------------------------------------------------
#> env_indices -- builds interval-resolved environmental indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ------------------------------------------------------------
#> env_indices(): 5 environments x 21 factors x 44 intervals = 924 indices
#>   interval width : 3 day(s), days 1..130
#>   discarded      : 106 record(s) outside 1..130
#> ------------------------------------------------------------

## Weekly windows instead
W3 <- env_indices(maizeWTH, interval = 7, end.day = 140)
#> ---------------------------------------------------------------
#> env_indices -- builds interval-resolved environmental indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ------------------------------------------------------------
#> env_indices(): 5 environments x 21 factors x 20 intervals = 420 indices
#>   interval width : 7 day(s), days 1..140
#>   discarded      : 60 record(s) outside 1..140
#> ------------------------------------------------------------
# }
```
