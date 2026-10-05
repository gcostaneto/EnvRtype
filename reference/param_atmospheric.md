# Estimate Atmospheric Parameters Related to Evapotranspiration

Computes vapour pressure deficit, the slope of the saturation vapour
pressure curve, Priestley-Taylor potential evapotranspiration and the
precipitation balance, following FAO-56.

## Usage

``` r
param_atmospheric(
  env.data,
  PREC = NULL,
  Tdew = NULL,
  Tmin = NULL,
  Tmax = NULL,
  RH = NULL,
  Rad = NULL,
  G = 0,
  Alt = 600,
  alpha = 1.26,
  merge = FALSE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. A `get_weather` output or equivalent.

- PREC, Tdew, Tmin, Tmax, RH, Rad:

  character. Column names for precipitation, dew point, minimum and
  maximum temperature, relative humidity and solar radiation. Defaults
  follow NASA POWER naming.

- G:

  numeric. Soil heat flux (MJ/m2/day). Default 0.

- Alt:

  numeric scalar, or character naming an elevation column in `env.data`.
  Default 600 m.

- alpha:

  numeric. Priestley-Taylor coefficient. Default 1.26.

- merge:

  boolean. Bind results onto `env.data`.

- verbose:

  boolean. Print a summary of what was computed.

## Value

A data.frame with `VPD` (kPa), `SPV` (kPa/C), `ETP` (mm/day) and `PETP`
(mm/day).

## Details

**Elevation now works.** `Alt` accepts either a numeric elevation or the
name of a column. Previously the default `Alt = 600` was never `NULL`,
so the elevation column was never read and every site was treated as 600
m – a 26% error in atmospheric pressure at 2500 m, propagating into the
psychrometric constant and ETP.

`VPD` is the difference between mean saturation vapour pressure (from
Tmin and Tmax) and actual vapour pressure (from dew point). Values are
floored at 0; negatives can only arise from a dew point above air
temperature, which is unphysical and indicates a data problem. A warning
is raised when this occurs.

## See also

`param_radiation`, `param_temperature`, `processWTH`

## Author

Germano Costa Neto

## Examples

``` r
if (FALSE) { # \dontrun{
env.data <- get_weather(lat = -13.05, lon = -56.05)
param_atmospheric(env.data, Alt = 850)
} # }
```
