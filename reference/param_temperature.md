# Estimate Thermal Parameters

Computes growing degree days, the temperature-limited radiation use
efficiency factor, and daily temperature range.

## Usage

``` r
param_temperature(
  env.data,
  Tmax = NULL,
  Tmin = NULL,
  Tbase1 = 9,
  Tbase2 = 45,
  Topt1 = 26,
  Topt2 = 32,
  method = c("standard", "legacy"),
  merge = FALSE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. A `get_weather` output or equivalent.

- Tmax, Tmin:

  character. Column names for maximum and minimum air temperature
  (Celsius).

- Tbase1:

  numeric. Base temperature below which development stops.

- Tbase2:

  numeric. Upper temperature at which development stops.

- Topt1, Topt2:

  numeric. Lower and upper bounds of the optimal range.

- method:

  character. GDD formulation: `"standard"` (default) or `"legacy"`. See
  Details.

- merge:

  boolean. Bind results onto `env.data`.

- verbose:

  boolean. Print a summary of what was computed.

## Value

A data.frame with `GDD` (C/day), `FRUE` (0-1) and `T2M_RANGE` (C).

## Details

**GDD methods.**

- `"standard"`:

  \\GDD = \max(\frac{\min(T\_{max}, T\_{base2}) + T\_{min}}{2} -
  T\_{base1},\\ 0)\\. Tmax is capped at the upper threshold; the result
  is floored at zero. This is the conventional agronomic definition.

- `"legacy"`:

  Clamps *both* Tmax and Tmin into \\\[T\_{base1}, T\_{base2}\]\\ before
  averaging. This is the behaviour of earlier EnvRtype versions and
  corresponds to the "Method 2" variant. It credits heat units on days
  colder than \\T\_{base1}\\ – e.g. with Tmax 12, Tmin 5, Tbase1 9 it
  returns 1.5 where the standard method returns 0.

The default changed from legacy to standard; results will differ on cold
days. Pass `method = "legacy"` to reproduce earlier output.

`FRUE` rises linearly from 0 at `Tbase1` to 1 at `Topt1`, holds at 1
across the optimum, and falls to 0 at `Tbase2`.

## See also

`param_radiation`, `param_atmospheric`, `processWTH`

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
env.data <- get_weather(lat = -13.05, lon = -56.05)
#> ---------------------------------------------------------------
#> get_weather -- downloads daily weather from NASA POWER
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> start.day is NULL; matched as 2024-01-18
#> end.day is NULL; matched as 2024-02-17
#> ---------------------------------------------------------------
#> get_weather() - daily weather from NASA POWER
#> https://docs.ropensci.org/nasapower
#> ---------------------------------------------------------------
#> Start .............................. Wed Oct 14 03:16:26 2026
#> Environmental units ................ 1
#> Parallel ........................... [   ]
#>   [1/1] env1 downloaded
#> Note: no valid data returned for: FROST_DAYS
#>       These columns are entirely NA (parameter unavailable at this site/period).
#> Computing VPD ...................... [ x ]
#> Computing N, RTA ................... [ x ]
#> Computing n ........................ [ x ]
#> Computing TH1 ...................... [ x ]
#> Computing TH2 ...................... [ x ]
#> Computing PAR_TEMP ................. [ x ]
#> Environmental units downloaded ..... 1/1
#> Daily weather features ............. 25
#> Done ............................... Wed Oct 14 03:16:27 2026
#> Total time ......................... 0.92s
param_temperature(env.data)
#> ---------------------------------------------------------------
#> param_temperature -- derives temperature-response parameters
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ----------------------------------------------------------------------
#> Growing Degree Day (GDD, C/day) - method: standard
#> Effect of temperature on radiation use efficiency (FRUE, 0-1)
#> Daily temperature range (T2M_RANGE, C)
#> ----------------------------------------------------------------------
#>       GDD      FRUE T2M_RANGE
#> 1  18.990 1.0000000     10.22
#> 2  17.955 1.0000000      8.15
#> 3  17.810 1.0000000      9.14
#> 4  17.825 1.0000000      6.95
#> 5  19.125 1.0000000      9.67
#> 6  18.495 1.0000000      7.43
#> 7  18.195 1.0000000      7.05
#> 8  18.360 1.0000000      6.58
#> 9  18.335 1.0000000     10.03
#> 10 18.225 1.0000000     11.21
#> 11 18.430 1.0000000     11.12
#> 12 17.790 1.0000000      9.00
#> 13 18.670 1.0000000     10.94
#> 14 18.455 1.0000000      9.05
#> 15 19.595 1.0000000     12.83
#> 16 17.355 1.0000000      6.11
#> 17 16.715 0.9832353      6.73
#> 18 15.595 0.9173529      5.21
#> 19 14.975 0.8808824      5.41
#> 20 16.105 0.9473529      8.17
#> 21 17.415 1.0000000      9.45
#> 22 17.370 1.0000000      9.38
#> 23 17.635 1.0000000      8.49
#> 24 17.805 1.0000000      8.05
#> 25 16.505 0.9708824      4.45
#> 26 16.505 0.9708824      6.23
#> 27 16.880 0.9929412      9.44
#> 28 17.705 1.0000000      9.69
#> 29 18.060 1.0000000      9.92
#> 30 18.095 1.0000000     10.17
#> 31 18.510 1.0000000      9.88
param_temperature(env.data, method = "legacy")  # pre-refactor behaviour
#> ---------------------------------------------------------------
#> param_temperature -- derives temperature-response parameters
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ----------------------------------------------------------------------
#> Growing Degree Day (GDD, C/day) - method: legacy
#> Effect of temperature on radiation use efficiency (FRUE, 0-1)
#> Daily temperature range (T2M_RANGE, C)
#> ----------------------------------------------------------------------
#>       GDD      FRUE T2M_RANGE
#> 1  18.990 1.0000000     10.22
#> 2  17.955 1.0000000      8.15
#> 3  17.810 1.0000000      9.14
#> 4  17.825 1.0000000      6.95
#> 5  19.125 1.0000000      9.67
#> 6  18.495 1.0000000      7.43
#> 7  18.195 1.0000000      7.05
#> 8  18.360 1.0000000      6.58
#> 9  18.335 1.0000000     10.03
#> 10 18.225 1.0000000     11.21
#> 11 18.430 1.0000000     11.12
#> 12 17.790 1.0000000      9.00
#> 13 18.670 1.0000000     10.94
#> 14 18.455 1.0000000      9.05
#> 15 19.595 1.0000000     12.83
#> 16 17.355 1.0000000      6.11
#> 17 16.715 0.9832353      6.73
#> 18 15.595 0.9173529      5.21
#> 19 14.975 0.8808824      5.41
#> 20 16.105 0.9473529      8.17
#> 21 17.415 1.0000000      9.45
#> 22 17.370 1.0000000      9.38
#> 23 17.635 1.0000000      8.49
#> 24 17.805 1.0000000      8.05
#> 25 16.505 0.9708824      4.45
#> 26 16.505 0.9708824      6.23
#> 27 16.880 0.9929412      9.44
#> 28 17.705 1.0000000      9.69
#> 29 18.060 1.0000000      9.92
#> 30 18.095 1.0000000     10.17
#> 31 18.510 1.0000000      9.88
# }
```
