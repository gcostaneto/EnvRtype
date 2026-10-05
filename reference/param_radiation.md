# Estimate Basic Solar-Radiation Parameters

Computes extraterrestrial radiation, daylength and actual sunshine
duration from day of year and latitude, following FAO-56.

## Usage

``` r
param_radiation(
  env.data,
  day.id = "DOY",
  latitude = "LAT",
  solarRadiation = NULL,
  merge = TRUE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. A `get_weather` output or equivalent.

- day.id:

  character. Column holding day of year. Default `"DOY"`.

- latitude:

  character. Column holding latitude. Default `"LAT"`.

- solarRadiation:

  character. Column holding incoming shortwave radiation. Default
  `"ALLSKY_SFC_SW_DWN"`.

- merge:

  boolean. Bind results onto `env.data`. Default `TRUE`.

- verbose:

  boolean. Print a summary of what was computed.

## Value

A data.frame with `n` (actual sunshine hours), `N` (daylength, h) and
`RTA` (extraterrestrial radiation, MJ/m2/day).

## Details

Daylength uses Forsythe et al. (1995), consistent with `get_weather`.
The two functions previously used different formulas that disagreed by
up to 0.5 h at high latitude.

The sunset hour angle is clamped so polar day/night yields 0 rather than
`NaN`. Missing radiation is median-imputed within the supplied data, and
`n` is bounded to `[0, N]`.

## See also

`param_atmospheric`, `param_temperature`, `processWTH`

## Author

Germano Costa Neto

## Examples

``` r
if (FALSE) { # \dontrun{
env.data <- get_weather(lat = -13.05, lon = -56.05)
param_radiation(env.data)
} # }
```
