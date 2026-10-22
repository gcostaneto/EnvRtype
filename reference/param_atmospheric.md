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
# \donttest{
env.data <- get_weather(lat = -13.05, lon = -56.05)
#> ---------------------------------------------------------------
#> get_weather -- downloads daily weather from NASA POWER
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> start.day is NULL; matched as 2024-01-26
#> end.day is NULL; matched as 2024-02-25
#> ---------------------------------------------------------------
#> get_weather() - daily weather from NASA POWER
#> https://docs.ropensci.org/nasapower
#> ---------------------------------------------------------------
#> Start .............................. Thu Oct 22 18:34:50 2026
#> Environmental units ................ 1
#> Parallel ........................... [   ]
#>   [1/1] env1 downloaded
#> Computing VPD ...................... [ x ]
#> Computing N, RTA ................... [ x ]
#> Computing n ........................ [ x ]
#> Computing TH1 ...................... [ x ]
#> Computing TH2 ...................... [ x ]
#> Computing PAR_TEMP ................. [ x ]
#> Environmental units downloaded ..... 1/1
#> Daily weather features ............. 25
#> Done ............................... Thu Oct 22 18:34:51 2026
#> Total time ......................... 0.96s
param_atmospheric(env.data, Alt = 850)
#> ---------------------------------------------------------------
#> param_atmospheric -- derives atmospheric-demand parameters
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ----------------------------------------------------------------------
#> Slope of saturation vapour pressure curve (SPV, kPa/C)
#> Vapour pressure deficit (VPD, kPa)
#> Potential evapotranspiration (ETP, mm/day)
#> Precipitation balance P - ETP (PETP, mm/day)
#> ----------------------------------------------------------------------
#>          VPD       SPV       ETP        PETP
#> 1  1.0570758 0.2127664 10.701635  -9.7216352
#> 2  1.1963363 0.2115765 11.262974 -11.2229737
#> 3  1.2402607 0.2137985 10.129089 -10.1090887
#> 4  0.9479085 0.2069256  8.669066  -4.7790658
#> 5  1.1868980 0.2164249  9.598868   1.5911318
#> 6  0.9927726 0.2140709  8.947486  -1.8774863
#> 7  1.3280358 0.2268024  9.526170  11.1938296
#> 8  0.6186276 0.2023610  6.961476   1.2085238
#> 9  0.7202045 0.1957995  4.285090   6.0149101
#> 10 0.3776662 0.1847477  5.813851  25.0261486
#> 11 0.3128784 0.1788592  7.653222   4.7067785
#> 12 0.6340613 0.1897132  9.404194   2.6958056
#> 13 0.7824328 0.2029855  8.081383  -3.5813835
#> 14 0.8544113 0.2025170  8.436664   5.5833361
#> 15 0.6804262 0.2052893  9.529413   1.0205868
#> 16 0.6819151 0.2070845  5.632104  13.2778964
#> 17 0.3567582 0.1936859  6.436445  15.8835547
#> 18 0.5607959 0.1936859  5.106924  -4.1769244
#> 19 0.6908651 0.1974738  9.106010  -7.6160103
#> 20 0.7864117 0.2060269  9.878412  -9.1184117
#> 21 0.8475449 0.2098021 10.333451  -5.6034514
#> 22 0.8393412 0.2101774  9.779679  -4.2996790
#> 23 0.9634764 0.2146710  6.338395  -3.3683950
#> 24 0.4906680 0.1902060  7.794400  22.3855999
#> 25 0.3687032 0.1950424  5.491460  15.0385401
#> 26 0.5029901 0.1938363  6.300752  10.2092476
#> 27 0.5504202 0.1993139  6.708586  -3.1685857
#> 28 0.6103149 0.2042394  5.962741  -4.2327406
#> 29 0.5414807 0.1995707  4.985698  -3.4756975
#> 30 0.5193414 0.1956984  7.502888  -0.7828879
#> 31 0.7193274 0.2021532  9.547041  -7.2970405
# }
```
