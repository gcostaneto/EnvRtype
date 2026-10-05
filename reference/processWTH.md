# Enrich a get_weather() Output with Derived Parameters

Convenience wrapper chaining `param_radiation`, `param_atmospheric` and
`param_temperature`.

## Usage

``` r
processWTH(
  env.data,
  PREC = NULL,
  Tdew = NULL,
  Tmin = NULL,
  Tmax = NULL,
  RH = NULL,
  Rad = NULL,
  G = 0,
  Alt = 600,
  Tbase1 = 9,
  Tbase2 = 45,
  Topt1 = 26,
  Topt2 = 32,
  method = c("standard", "legacy"),
  alpha = 1.26,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. A `get_weather` output.

- PREC, Tdew, Tmin, Tmax, RH, Rad:

  character. Column names passed through to the underlying functions.

- G:

  numeric. Soil heat flux. Default 0.

- Alt:

  numeric or character. Elevation in metres, or the name of an elevation
  column.

- Tbase1, Tbase2, Topt1, Topt2:

  numeric. Cardinal temperatures.

- method:

  character. GDD method, see `param_temperature`.

- alpha:

  numeric. Priestley-Taylor coefficient.

- verbose:

  boolean. Print progress messages.

## Value

The input data.frame with radiation, atmospheric and thermal parameters
appended.

## Details

If `env.data` already carries `n`, `N` or `RTA` from `get_weather`,
those columns are dropped before recomputation to avoid duplicated
names. Note that
[`get_weather()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md)
derives `n` from `ALLSKY_SFC_SW_DNI` whereas
[`param_radiation()`](https://gcostaneto.github.io/EnvRtype/reference/param_radiation.md)
uses `ALLSKY_SFC_SW_DWN`, so the recomputed values may differ.

## See also

`get_weather`, `summaryWTH`

## Author

Germano Costa Neto

## Examples

``` r
if (FALSE) { # \dontrun{
env.data <- get_weather(lat = -13.05, lon = -56.05)
processWTH(env.data)
} # }
```
