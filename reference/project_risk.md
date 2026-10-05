# Compare Envirotype Risk Between Baseline and a Climate Scenario

Runs
[`env_risk_profile()`](https://gcostaneto.github.io/EnvRtype/reference/env_risk_profile.md)
on an observed series and on a delta-shifted scenario series, and
returns the change in stress frequency per site, stage and class.

## Usage

``` r
project_risk(baseline, scenario, ..., verbose = TRUE)
```

## Arguments

- baseline:

  data.frame. Observed multi-year daily weather.

- scenario:

  data.frame. Output of
  [`get_climate_scenario()`](https://gcostaneto.github.io/EnvRtype/reference/get_climate_scenario.md)
  built from the same baseline, or any comparably structured series.

- ...:

  Arguments passed to
  [`env_risk_profile()`](https://gcostaneto.github.io/EnvRtype/reference/env_risk_profile.md)
  for BOTH runs, so the two profiles are strictly comparable.

- verbose:

  boolean. Print progress.

## Value

A list of class `"env_risk_delta"` with `baseline`, `scenario` (both
`"env_risk"` objects) and `change`, a data.frame with `freq_base`,
`freq_scen`, `freq_diff` and `days_diff`.

## Details

Because `get_climate_scenario(method = "delta")` preserves the observed
day-to-day sequence, differences here reflect the MEAN shift only.
Rising frequencies of heat classes are real but conservative; unchanged
precipitation-class frequencies do not mean rainfall risk is unchanged,
only that the mean-preserving perturbation left wet-day counts alone.

## See also

[`env_risk_profile`](https://gcostaneto.github.io/EnvRtype/reference/env_risk_profile.md),
[`get_climate_scenario`](https://gcostaneto.github.io/EnvRtype/reference/get_climate_scenario.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fut <- get_climate_scenario(env.id = sites$env, lat = sites$lat,
                            lon = sites$lon, baseline = wth,
                            scenario = "ssp585", period = "2061-2080")

chg <- project_risk(baseline = processWTH(wth),
                    scenario = processWTH(fut),
                    planting = "10-15", crop = "maize",
                    var.id = c("T2M_MAX", "PETP"))
subset(chg$change, class == "heat_severe")
} # }
```
