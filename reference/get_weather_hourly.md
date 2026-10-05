# Day/Night Partitioned Hourly Weather from NASA POWER

Sibling collector to
[`get_weather`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md)
that downloads the NASA POWER *hourly* endpoint, partitions each 24 h
cycle into day / night / pre-dawn, and returns one row per
environment-day. It is a PARALLEL producer of a covariate block; it does
NOT feed
[`processWTH()`](https://gcostaneto.github.io/EnvRtype/reference/processWTH.md)
and must not be wired into the daily pipeline by analogy.

## Usage

``` r
get_weather_hourly(
  env.id,
  lat,
  lon,
  start.day,
  end.day,
  variables.names = NULL,
  night.def = c("solar", "fixed", "radiation"),
  night.hours = NULL,
  predawn.hours = 2:5,
  time_standard = c("LST", "UTC"),
  summarise = TRUE,
  thresholds = NULL,
  tries = 3,
  save = FALSE,
  dir.path = NULL,
  verbose = TRUE
)
```

## Arguments

- env.id:

  character. Environment identifiers.

- lat, lon:

  numeric. Coordinates, recycled if length 1.

- start.day, end.day:

  character/Date. Collection window, recycled if length 1.

- variables.names:

  character. POWER parameters. **Maximum 15** on the hourly endpoint;
  validated before the request is sent.

- night.def:

  character. How to define night:

  `"solar"`

  :   Sunset/sunrise from Forsythe daylength. Night length varies with
      season and latitude.

  `"fixed"`

  :   User window via `night.hours`. Constant and directly comparable
      across sites and dates.

  `"radiation"`

  :   Hours with `ALLSKY_SFC_SW_DWN <= 0`. Empirical; costs one
      parameter slot.

- night.hours:

  integer. Hours 0-23 for `night.def = "fixed"`. Wrap-around (e.g.
  `c(22,23,0,1,2,3,4,5)`) is handled.

- predawn.hours:

  integer. Window for pre-dawn statistics. Default 2:5.

- time_standard:

  character. `"LST"` (default) or `"UTC"`. Keep LST: hour 3 is then ~3
  a.m. solar time at every site, so pre-dawn windows are comparable
  across longitudes.

- summarise:

  boolean. `TRUE` (default) returns one row per environment-day. `FALSE`
  returns raw hourly rows plus `period` and `night_id` – the 24x volume
  path.

- thresholds:

  named numeric. Night-hour counts above each threshold, e.g.
  `c(T2M = 22)` gives `hours_above_T2M`.

- tries:

  integer. Retry attempts per site.

- save:

  boolean. Write one per-site CSV.

- dir.path:

  character. Output directory when `save = TRUE`.

- verbose:

  boolean. Print progress messages.

## Value

A data.frame of class `c("weather_hourly", "data.frame")`. With
`summarise = TRUE`, one row per environment-day, carrying day/night/
pre-dawn statistics per variable. **Radiation-type parameters
(shortwave, longwave, PAR) are summarised for DAYTIME only**: they are
~0 at night by construction, so night statistics would add no
information.

## Details

**POWER hourly is MERRA-2 reanalysis, not observation.** The diurnal
cycle at 0.5 x 0.625 degrees is modelled and sub-daily fidelity is
weaker than the daily aggregate.

A daily `T2M_MIN` is not a substitute for a night mean: the two differ
materially, and on a `Q10 ~ 2` process a few degrees is a large
difference in implied rate (Peng et al. 2004).

## References

Peng, S. et al. (2004) Rice yields decline with higher night temperature
from global warming. *PNAS* 101(27), 9971-9975.

Forsythe, W.C. et al. (1995) A model comparison for daylength as a
function of latitude and day of the year. *Ecological Modelling* 80,
87-95.

Allen, R.G. et al. (1998) Crop evapotranspiration. *FAO Irrigation and
Drainage Paper* 56.

## See also

[`get_weather`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md),
[`processWTH`](https://gcostaneto.github.io/EnvRtype/reference/processWTH.md),
[`summaryWTH`](https://gcostaneto.github.io/EnvRtype/reference/summaryWTH.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Reuse the coordinates of a shipped daily dataset
sites <- unique(maizeWTH[, c("env", "LON", "LAT")])
h <- get_weather_hourly(env.id = sites$env, lat = sites$LAT, lon = sites$LON,
                        start.day = "2016-01-01", end.day = "2016-01-31",
                        night.def = "solar")

## Night temperature is not the daily minimum
with(h, cor(T2M_night_mean, T2M_day_min, use = "complete.obs"))

## Fixed night window, wrap-around handled correctly
h2 <- get_weather_hourly(sites$env, sites$LAT, sites$LON,
                         "2016-01-01", "2016-01-31",
                         night.def = "fixed", night.hours = c(22:23, 0:5))
} # }
```
