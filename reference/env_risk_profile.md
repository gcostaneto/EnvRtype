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
# \donttest{
## 30 seasons of weather for three sites
sites <- data.frame(env = c("PIRA", "SETE", "NMAC"),
                    lat = c(-22.7, -19.4, -21.2),
                    lon = c(-47.6, -44.2, -45.2))

wth <- get_weather(env.id = sites$env, lat = sites$lat, lon = sites$lon,
                   start.day = "1994-01-01", end.day = "2024-12-31")
#> ---------------------------------------------------------------
#> get_weather -- downloads daily weather from NASA POWER
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ---------------------------------------------------------------
#> get_weather() - daily weather from NASA POWER
#> https://docs.ropensci.org/nasapower
#> ---------------------------------------------------------------
#> Start .............................. Thu Oct 22 18:22:51 2026
#> Environmental units ................ 3
#> Parallel ........................... [   ]
#>   [1/3] PIRA downloaded
#>   [2/3] SETE downloaded
#>   [3/3] NMAC downloaded
#> Computing VPD ...................... [ x ]
#> Computing N, RTA ................... [ x ]
#> Computing n ........................ [ x ]
#> Computing TH1 ...................... [ x ]
#> Computing TH2 ...................... [ x ]
#> Computing PAR_TEMP ................. [ x ]
#> Environmental units downloaded ..... 3/3
#> Daily weather features ............. 25
#> Done ............................... Thu Oct 22 18:22:58 2026
#> Total time ......................... 6.57s
wth <- processWTH(wth)
#> ---------------------------------------------------------------
#> processWTH -- derives all agro-meteorological parameters
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Dropping existing column(s) before recomputation: n, N, RTA, VPD
#>   - radiation parameters computed
#>   - atmospheric-demand parameters computed
#>   - temperature-response parameters computed

## Risk of heat and drought at flowering, sown 15 October each year
rp <- env_risk_profile(wth, site.id = "env", planting = "10-15",
                       var.id = c("T2M_MAX", "PETP"),
                       crop = "maize", by.stage = TRUE)
#> ---------------------------------------------------------------
#> env_risk_profile -- multi-year envirotype frequency and risk
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - assembling 3 site(s) x 31 season(s)
#>   - 93 complete site-year(s) retained
#>   - assigning thermal-time stages per season
#> Warning: 93 environment(s) did not accumulate enough thermal time to reach 'R6': NMAC_1994, NMAC_1995, NMAC_1996, NMAC_1997, NMAC_1998, NMAC_1999, NMAC_2000, NMAC_2001, .... Days beyond the last reached stage are NA. Extend the series or lower the stage thresholds.
#>   - classifying stress envirotypes
rp
#> --------------------------------------------------------------- 
#> env_risk_profile() - historical envirotype risk
#> --------------------------------------------------------------- 
#> Sites .............. 3
#> Seasons ............ 31 (1994-2024)
#> Site-years ......... 93
#> Variables .......... PETP, T2M_MAX
#> Stage-resolved ..... [ x ]
#> 
#> Most frequent stress classes:
#>   SETE     R3     PETP      drought_severe freq 1.00   16.4 d/season
#>   PIRA     V6     PETP      drought_severe freq 1.00   15.6 d/season
#>   NMAC     R3     PETP      drought_severe freq 1.00   15.3 d/season
#>   PIRA     R5     PETP      drought_severe freq 1.00   14.5 d/season
#>   NMAC     V6     PETP      drought_severe freq 1.00   14.4 d/season
#>   SETE     R5     PETP      drought_severe freq 1.00   14.1 d/season
#>   SETE     V6     PETP      drought_severe freq 1.00   13.7 d/season
#>   PIRA     R3     PETP      drought_severe freq 1.00   13.2 d/season
#>   PIRA     VT     PETP      drought_severe freq 1.00   12.9 d/season
#>   PIRA     V6     T2M_MAX   heat_mild      freq 1.00   12.5 d/season
#> --------------------------------------------------------------- 
#> $frequency  $severity  $return_period

## How often is flowering heat-stressed at each site?
subset(rp$frequency, stage == "R1" & class == "heat_severe")
#>    site stage variable       class n_years n_occur       freq  mean_days
#> 2  NMAC    R1  T2M_MAX heat_severe      30       1 0.03333333 0.23333333
#> 23 PIRA    R1  T2M_MAX heat_severe      30       1 0.03333333 0.03333333
#> 47 SETE    R1  T2M_MAX heat_severe      30       2 0.06666667 0.13333333
#>    max_days
#> 2         7
#> 23        1
#> 47        2

## Weights for a target population of environments
w <- tpe_weights(rp, stage = "R1", var = "PETP")
#> ------------------------------------------------ 
#> TPE weights by class (measure = freq)
#> ------------------------------------------------ 
#>   drought_severe      0.364  ===========
#>   drought_mild        0.329  ==========
#>   wet                 0.307  =========
#> ------------------------------------------------ 
# }
```
