# Basic Summary Statistics for Environmental Data

Summarises a `get_weather` output by environment and, optionally, by
user-defined time intervals such as phenological stages.

## Usage

``` r
summaryWTH(
  env.data,
  id.names = NULL,
  env.id = NULL,
  days.id = NULL,
  var.id = NULL,
  statistic = NULL,
  probs = NULL,
  by.interval = FALSE,
  time.window = NULL,
  names.window = NULL,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. A `get_weather` output.

- id.names:

  vector (character). Columns treated as identifiers.

- env.id:

  character. Environment-id column. Default `"env"`.

- days.id:

  character. Column holding days from start. Default `"daysFromStart"`.

- var.id:

  vector (character). Variables to summarise. Defaults to every
  non-identifier column.

- statistic:

  character. One of `"all"`, `"sum"`, `"mean"`, `"quantile"`.

- probs:

  vector (numeric). Quantile probabilities. Default
  `c(0.25, 0.50, 0.75)`.

- by.interval:

  boolean. Summarise within time intervals.

- time.window:

  vector (numeric). Interval breakpoints in days.

- names.window:

  vector (character). Interval names.

- verbose:

  boolean. If `TRUE` (default) prints a progress banner.

## Value

A data.frame in long format with one row per environment (and interval)
and variable.

## Details

Quantiles are computed with
[`stats::aggregate`](https://rdrr.io/r/stats/aggregate.html) rather than
the previous nested `foreach` loops, which registered no parallel
backend and therefore ran sequentially while emitting a warning on every
call.

## See also

`get_weather`, `processWTH`

## Author

Germano Costa Neto

## Examples

``` r
if (FALSE) { # \dontrun{
env.data <- get_weather(lat = -13.05, lon = -56.05)
summaryWTH(env.data, env.id = "env", statistic = "mean")

summaryWTH(env.data, env.id = "env", by.interval = TRUE,
           time.window  = c(0, 14, 35, 60, 90, 120),
           names.window = c("P-E", "E-V1", "V1-V4", "V4-VT", "VT-GF", "GF-PM"))
} # }
```
