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
#> start.day is NULL; matched as 2024-01-18
#> end.day is NULL; matched as 2024-02-17
#> ---------------------------------------------------------------
#> get_weather() - daily weather from NASA POWER
#> https://docs.ropensci.org/nasapower
#> ---------------------------------------------------------------
#> Start .............................. Wed Oct 14 03:16:24 2026
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
#> Done ............................... Wed Oct 14 03:16:25 2026
#> Total time ......................... 1.17s
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
#>          VPD       SPV       ETP       PETP
#> 1  1.1329342 0.2199688 10.679594  -6.419594
#> 2  0.8908138 0.2086795  8.653845  19.106155
#> 3  0.7845527 0.2071375 10.422778  -2.872778
#> 4  0.6764902 0.2072966  9.653876  -7.783876
#> 5  1.0276807 0.2214785  7.748428  -4.998428
#> 6  0.8150283 0.2145072  7.069940  -2.239940
#> 7  0.6079081 0.2112529  7.899735   0.580265
#> 8  0.7429793 0.2130376  6.959234   2.210766
#> 9  1.0570758 0.2127664 10.701635  -9.721635
#> 10 1.1963363 0.2115765 11.262974 -11.222974
#> 11 1.2402607 0.2137985 10.129089 -10.109089
#> 12 0.9479085 0.2069256  8.669066  -4.779066
#> 13 1.1868980 0.2164249  9.598868   1.591132
#> 14 0.9927726 0.2140709  8.947486  -1.877486
#> 15 1.3280358 0.2268024  9.526170  11.193830
#> 16 0.6186276 0.2023610  6.961476   1.208524
#> 17 0.7202045 0.1957995  4.285090   6.014910
#> 18 0.3776662 0.1847477  5.813851  25.026149
#> 19 0.3128784 0.1788592  7.653222   4.706778
#> 20 0.6340613 0.1897132  9.404194   2.695806
#> 21 0.7824328 0.2029855  8.081383  -3.581383
#> 22 0.8544113 0.2025170  8.436664   5.583336
#> 23 0.6804262 0.2052893  9.529413   1.020587
#> 24 0.6819151 0.2070845  5.632104  13.277896
#> 25 0.3567582 0.1936859  6.436445  15.883555
#> 26 0.5607959 0.1936859  5.106924  -4.176924
#> 27 0.6908651 0.1974738  9.106010  -7.616010
#> 28 0.7864117 0.2060269  9.878412  -9.118412
#> 29 0.8475449 0.2098021 10.333451  -5.603451
#> 30 0.8393412 0.2101774  9.779679  -4.299679
#> 31 0.9634764 0.2146710  6.338395  -3.368395
# }
```
