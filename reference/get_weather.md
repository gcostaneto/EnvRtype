# Easily Collect Worldwide Daily Weather Data

Imports daily weather data from the NASA POWER API and augments it with
derived agro-meteorological indices (VPD, photoperiod, extraterrestrial
radiation, temperature-humidity indices).

## Usage

``` r
get_weather(
  env.id = NULL,
  lat = NULL,
  lon = NULL,
  start.day = NULL,
  end.day = NULL,
  variables.names = NULL,
  dir.path = NULL,
  save = FALSE,
  parallel = FALSE,
  workers = NULL,
  chunk_size = 29,
  sleep = 60,
  tries = 3L,
  verbose = TRUE
)
```

## Arguments

- env.id:

  vector (character). Identifier of the site/environment. If `NULL`,
  environments are named `env1 ... envN`.

- lat:

  vector (numeric). Latitude in decimal degrees (WGS84).

- lon:

  vector (numeric). Longitude in decimal degrees (WGS84).

- start.day, end.day:

  vector (character or Date). First and last date of collection, e.g.
  `"2015-02-15"`. Recycled if length 1.

- variables.names:

  vector (character). POWER parameters to request. See Details.

- dir.path:

  character. Output directory when `save = TRUE`. Defaults to
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- save:

  boolean. Write one .csv per environment.

- parallel:

  boolean. Download chunks in parallel. Default `FALSE` (serial), which
  is gentler on the NASA POWER rate limit; set `TRUE` to speed up large
  jobs.

- workers:

  integer. Number of parallel processes. Defaults to 90% of available
  cores, capped at `chunk_size`.

- chunk_size:

  integer. Points per chunk. Default 29; raising this may exceed the API
  rate limit.

- sleep:

  numeric. Seconds to pause between chunks. Default 60.

- tries:

  integer. Retry attempts per failed request. Default 3.

- verbose:

  boolean. Print progress messages.

## Value

A data.frame with one row per environment-day. Always contains `env`,
`LON`, `LAT`, `YYYYMMDD`, `DOY` and `daysFromStart`, followed by the
requested and derived variables.

## Details

**daysFromStart is date-based.** It is computed as
`as.numeric(YYYYMMDD - start.day) + 1`, not as the row position. The
previous implementation used `1:nrow()`, which silently mis-numbered
every row after a gap whenever the API returned fewer days than
requested – assigning observations to the wrong phenological stage in
`summaryWTH`.

**Missing values.** POWER fill values (-999 and -99) are converted to
`NA` before any derived variable is computed.

Commonly used parameters include `T2M`, `T2M_MAX`, `T2M_MIN`, `T2MDEW`,
`PRECTOTCORR` (returned as `PRECTOT`), `RH2M`, `WS2M`,
`ALLSKY_SFC_SW_DWN`, `ALLSKY_SFC_SW_DNI`, `ALLSKY_SFC_PAR_TOT`,
`GWETROOT`, `GWETTOP`, `EVPTRNS`.

Derived variables, added when their inputs are present:

- `VPD`: vapour pressure deficit (kPa), from T2MDEW, T2M_MAX, T2M_MIN

- `N`: photoperiod (h), Forsythe et al. (1995)

- `RTA`: extraterrestrial radiation (MJ/m2/day), FAO-56

- `n`: actual sunshine duration (h)

- `TH1`, `TH2`: temperature-humidity indices (NRC 1971; Yousef 1985)

- `PAR_TEMP`: ratio of PAR to mean temperature

## References

Sparks, A. (2018). nasapower: NASA-POWER Data from R. *JOSS* 3(30),
1035.

Forsythe, W.C., Rykiel, E.J., Stahl, R.S., Wu, H., Schoolfield, R.M.
(1995). A model comparison for daylength as a function of latitude and
day of the year. *Ecological Modelling* 80, 87-95.

## See also

`summaryWTH`, `processWTH`, `get_soil`,
[`get_weather_hourly`](https://gcostaneto.github.io/EnvRtype/reference/get_weather_hourly.md)

## Author

Germano Costa Neto and Giovanni Galli, modified by Tiago Olivoto

## Examples

``` r
if (FALSE) { # \dontrun{
## Single location
get_weather(env.id = "NM", lat = -13.05, lon = -56.05,
            start.day = "2015-02-15", end.day = "2015-06-15",
            variables.names = c("T2M", "PRECTOTCORR"))

## Two locations with different planting dates
get_weather(env.id = c("NM", "SO"),
            lat = c(-13.05, -12.32), lon = c(-56.05, -55.42),
            start.day = c("2015-02-15", "2015-02-13"),
            end.day   = rep("2015-06-15", 2))
} # }
```
