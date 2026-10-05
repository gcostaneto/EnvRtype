# Historical Envirotype Risk Profile Across Years

Quantifies how OFTEN each stress envirotype occurs at a site, by
developmental stage, across many historical seasons. Where
[`env_typing()`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md)
describes the environments you actually observed, `env_risk_profile()`
estimates the underlying climatology those environments were drawn from
– the frequency, not the realisation.

This is the quantity a breeding programme needs to weight its testing
network: a mega-environment that appears in 1 season out of 20 should
not carry the same weight as one that appears in 12.

## Usage

``` r
env_risk_profile(
  env.data,
  site.id = "env",
  date.id = "YYYYMMDD",
  var.id = NULL,
  years = NULL,
  planting = NULL,
  season.length = 150,
  cardinals = NULL,
  by.stage = TRUE,
  crop = "maize",
  stages = NULL,
  Tbase = NULL,
  Tupper = NULL,
  method = "gdd",
  min.days = 30L,
  return.seasons = FALSE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. Multi-year daily weather for one or more sites, typically
  [`get_weather()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md)
  run over a long window, then
  [`processWTH()`](https://gcostaneto.github.io/EnvRtype/reference/processWTH.md).
  Must contain a site id, a date, and the variables in `var.id`.

- site.id:

  character. Site-id column. Default `"env"`. Note this is a SITE, not a
  site-year; the site-year is formed internally.

- date.id:

  character. Date column (`Date` or `YYYYMMDD`). Default `"YYYYMMDD"`.

- var.id:

  vector (character). Variables to profile. Defaults to those present
  among the names of `cardinals`.

- years:

  vector (numeric). Seasons to include. Defaults to every year present
  with a complete season.

- planting:

  character, `Date` vector or `"DOY"`. Planting date per season. A
  single `"MM-DD"` string applies the same sowing date every year; a
  named vector keyed by site applies a per-site date; a `Date` vector is
  used as given.

- season.length:

  numeric. Days after planting retained per season when phenology is not
  used. Default 150.

- cardinals:

  named list. Per-variable `list(breaks =, labels =)` defining stress
  classes. Defaults in `.RISK_CARDINALS` cover `T2M_MAX`, `T2M_MIN`,
  `T2M`, `PETP`, `PRECTOT`, `VPD`, `FRUE` and `Ks`.

- by.stage:

  boolean. If `TRUE` (default) compute risk within thermal-time stages
  by calling
  [`env_phenology()`](https://gcostaneto.github.io/EnvRtype/reference/env_phenology.md)
  per season.

- crop, stages, Tbase, Tupper, method:

  Passed to
  [`env_phenology()`](https://gcostaneto.github.io/EnvRtype/reference/env_phenology.md)
  when `by.stage = TRUE`.

- min.days:

  integer. Minimum days a site-year must contribute to be kept. Default
  30.

- return.seasons:

  boolean. If `TRUE`, also return the per-site-year table behind the
  frequencies.

- verbose:

  boolean. Print progress.

## Value

An object of class `"env_risk"`: a list with

- frequency:

  data.frame of site x stage x variable x class with `n_years`,
  `n_occur`, `freq` (probability a season shows that class at all), and
  `mean_days` (expected days per season)

- severity:

  data.frame of site x stage x variable with the mean, sd and 10th/90th
  percentile of the variable across seasons

- return_period:

  data.frame with `1 / freq`, the expected number of seasons between
  occurrences of each stress class

- seasons:

  per-site-year detail if `return.seasons = TRUE`

- call, years, cardinals:

  provenance

## Details

**Two frequencies are reported and they answer different questions.**
`freq` is the proportion of SEASONS in which the class occurred at least
once – the probability of encountering that stress. `mean_days` is the
expected NUMBER OF DAYS per season – the intensity. A site can have
`freq = 1.0` for mild heat while averaging only 3 days of it; that is a
very different target from one averaging 40 days.

**Return period** is simply `1 / freq` and inherits all the usual
caveats: it is an average recurrence interval estimated from a finite
record, not a schedule. With 20 years of data a return period above ~10
seasons is barely distinguishable from noise. The function will not
print a return period longer than the record.

**Stationarity.** Frequencies computed over 1990-2024 describe that
period. If the site is warming, they are already an underestimate of
future heat risk. Use
[`project_risk()`](https://gcostaneto.github.io/EnvRtype/reference/project_risk.md)
with a scenario to see the shift.

## References

Chapman, S.C. et al. (2000) Using crop simulation to generate genotype
by environment interaction effects for sorghum in water-limited
environments. *Aust. J. Agric. Res.* 51(2), 209-221.

Heinemann, A.B. et al. (2015) Drought impact on rainfed common bean
production areas in Brazil. *Agric. For. Meteorol.* 225, 1-12.

## See also

[`env_phenology`](https://gcostaneto.github.io/EnvRtype/reference/env_phenology.md),
[`env_typing`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md),
[`tpe_weights`](https://gcostaneto.github.io/EnvRtype/reference/tpe_weights.md),
[`project_risk`](https://gcostaneto.github.io/EnvRtype/reference/project_risk.md)

## Author

Germano Costa Neto

## Examples

``` r
if (FALSE) { # \dontrun{
## 30 seasons of weather for three sites
sites <- data.frame(env = c("PIRA", "SETE", "NMAC"),
                    lat = c(-22.7, -19.4, -21.2),
                    lon = c(-47.6, -44.2, -45.2))

wth <- get_weather(env.id = sites$env, lat = sites$lat, lon = sites$lon,
                   start.day = "1994-01-01", end.day = "2024-12-31")
wth <- processWTH(wth)

## Risk of heat and drought at flowering, sown 15 October each year
rp <- env_risk_profile(wth, site.id = "env", planting = "10-15",
                       var.id = c("T2M_MAX", "PETP"),
                       crop = "maize", by.stage = TRUE)
rp

## How often is flowering heat-stressed at each site?
subset(rp$frequency, stage == "R1" & class == "heat_severe")

## Weights for a target population of environments
w <- tpe_weights(rp, stage = "R1", var = "PETP")
} # }
```
