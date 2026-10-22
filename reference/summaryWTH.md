# Basic Summary Statistics for Environmental Data

Summarises a `get_weather` output by environment and, optionally, by
user-defined time intervals such as phenological stages.

## Usage

``` r
summaryWTH(
  env.data,
  id.names = NULL,
  env.id = NULL,
  days.id = NULL,
  var.id = NULL,
  statistic = NULL,
  probs = NULL,
  by.interval = FALSE,
  time.window = NULL,
  names.window = NULL,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. A `get_weather` output.

- id.names:

  vector (character). Columns treated as identifiers.

- env.id:

  character. Environment-id column. Default `"env"`.

- days.id:

  character. Column holding days from start. Default `"daysFromStart"`.

- var.id:

  vector (character). Variables to summarise. Defaults to every
  non-identifier column.

- statistic:

  character. One of `"all"`, `"sum"`, `"mean"`, `"quantile"`.

- probs:

  vector (numeric). Quantile probabilities. Default
  `c(0.25, 0.50, 0.75)`.

- by.interval:

  boolean. Summarise within time intervals.

- time.window:

  vector (numeric). Interval breakpoints in days.

- names.window:

  vector (character). Interval names.

- verbose:

  boolean. If `TRUE` (default) prints a progress banner.

## Value

A data.frame in long format with one row per environment (and interval)
and variable.

## Details

Quantiles are computed with
[`stats::aggregate`](https://rdrr.io/r/stats/aggregate.html) rather than
the previous nested `foreach` loops, which registered no parallel
backend and therefore ran sequentially while emitting a warning on every
call.

## See also

`get_weather`, `processWTH`

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
#> Start .............................. Thu Oct 22 18:35:09 2026
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
#> Done ............................... Thu Oct 22 18:35:11 2026
#> Total time ......................... 1.08s
summaryWTH(env.data, env.id = "env", statistic = "mean")
#> ---------------------------------------------------------------
#> summaryWTH -- summarises weather by environment and interval
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>     env           variable        mean
#> 1  env1                T2M 25.92709677
#> 2  env1            T2M_MAX 30.52935484
#> 3  env1            T2M_MIN 22.16516129
#> 4  env1             T2MDEW 22.89483871
#> 5  env1  ALLSKY_SFC_LW_DWN 36.08806452
#> 6  env1  ALLSKY_SFC_SW_DWN 20.05838710
#> 7  env1  ALLSKY_SFC_SW_DNI 11.89387097
#> 8  env1 ALLSKY_SFC_PAR_TOT  9.72709677
#> 9  env1     ALLSKY_SFC_UVA  1.34967742
#> 10 env1     ALLSKY_SFC_UVB  0.04290323
#> 11 env1            PRECTOT  9.26806452
#> 12 env1            EVPTRNS  9.04677419
#> 13 env1               QV2M 18.13903226
#> 14 env1               RH2M 84.88129032
#> 15 env1           GWETROOT  0.74548387
#> 16 env1            GWETTOP  0.79032258
#> 17 env1         FROST_DAYS  0.00000000
#> 18 env1               WS2M  0.48419355
#> 19 env1                VPD  0.74065646
#> 20 env1                  N 12.57587097
#> 21 env1                RTA 39.78167742
#> 22 env1                  n  3.76003226
#> 23 env1                TH1 76.91020373
#> 24 env1                TH2 75.36923871
#> 25 env1           PAR_TEMP  0.37401389

summaryWTH(env.data, env.id = "env", by.interval = TRUE,
           time.window  = c(0, 14, 35, 60, 90, 120),
           names.window = c("P-E", "E-V1", "V1-V4", "V4-VT", "VT-GF", "GF-PM"))
#> ---------------------------------------------------------------
#> summaryWTH -- summarises weather by environment and interval
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>     env interval           variable        mean         sum  prob_0.25
#> 1  env1     E-V1  ALLSKY_SFC_LW_DWN 36.30000000  653.400000 36.0125000
#> 2  env1     E-V1 ALLSKY_SFC_PAR_TOT  9.26500000  166.770000  7.6075000
#> 3  env1     E-V1  ALLSKY_SFC_SW_DNI 10.48444444  188.720000  4.5450000
#> 4  env1     E-V1  ALLSKY_SFC_SW_DWN 18.98333333  341.700000 15.2525000
#> 5  env1     E-V1     ALLSKY_SFC_UVA  1.28444444   23.120000  1.0775000
#> 6  env1     E-V1     ALLSKY_SFC_UVB  0.04055556    0.730000  0.0300000
#> 7  env1     E-V1            EVPTRNS  9.76777778  175.820000  7.7625000
#> 8  env1     E-V1         FROST_DAYS  0.00000000    0.000000  0.0000000
#> 9  env1     E-V1           GWETROOT  0.76500000   13.770000  0.7525000
#> 10 env1     E-V1            GWETTOP  0.80833333   14.550000  0.7800000
#> 11 env1     E-V1                  N 12.50661111  225.119000 12.4587500
#> 12 env1     E-V1           PAR_TEMP  0.35741969    6.433554  0.2919944
#> 13 env1     E-V1            PRECTOT  9.17388889  165.130000  1.8600000
#> 14 env1     E-V1               QV2M 18.53722222  333.670000 18.3000000
#> 15 env1     E-V1               RH2M 86.65055556 1559.710000 83.3300000
#> 16 env1     E-V1                RTA 39.54622222  711.832000 39.3470000
#> 17 env1     E-V1                T2M 25.85833333  465.450000 25.1125000
#> 18 env1     E-V1             T2MDEW 23.25500000  418.590000 23.0800000
#> 19 env1     E-V1            T2M_MAX 30.22888889  544.120000 29.3300000
#> 20 env1     E-V1            T2M_MIN 22.22722222  400.090000 21.7625000
#> 21 env1     E-V1                TH1 77.01414081 1386.254535 76.3980104
#> 22 env1     E-V1                TH2 75.43013333 1357.742400 74.8568000
#> 23 env1     E-V1                VPD  0.64251066   11.565192  0.5248762
#> 24 env1     E-V1               WS2M  0.52722222    9.490000  0.4600000
#> 25 env1     E-V1                  n  3.31538889   59.677000  1.4360000
#> 26 env1      P-E  ALLSKY_SFC_LW_DWN 35.79461538  465.330000 35.3300000
#> 27 env1      P-E ALLSKY_SFC_PAR_TOT 10.36692308  134.770000  9.7600000
#> 28 env1      P-E  ALLSKY_SFC_SW_DNI 13.84538462  179.990000  6.6100000
#> 29 env1      P-E  ALLSKY_SFC_SW_DWN 21.54692308  280.110000 19.9600000
#> 30 env1      P-E     ALLSKY_SFC_UVA  1.44000000   18.720000  1.3600000
#> 31 env1      P-E     ALLSKY_SFC_UVB  0.04615385    0.600000  0.0400000
#> 32 env1      P-E            EVPTRNS  8.04846154  104.630000  5.9600000
#> 33 env1      P-E         FROST_DAYS  0.00000000    0.000000  0.0000000
#> 34 env1      P-E           GWETROOT  0.71846154    9.340000  0.6900000
#> 35 env1      P-E            GWETTOP  0.76538462    9.950000  0.7200000
#> 36 env1      P-E                  N 12.67176923  164.733000 12.6430000
#> 37 env1      P-E           PAR_TEMP  0.39699048    5.160876  0.3806796
#> 38 env1      P-E            PRECTOT  9.39846154  122.180000  3.8900000
#> 39 env1      P-E               QV2M 17.58769231  228.640000 17.1300000
#> 40 env1      P-E               RH2M 82.43153846 1071.610000 78.2500000
#> 41 env1      P-E                RTA 40.10769231  521.400000 40.0420000
#> 42 env1      P-E                T2M 26.02230769  338.290000 24.9200000
#> 43 env1      P-E             T2MDEW 22.39615385  291.150000 21.9900000
#> 44 env1      P-E            T2M_MAX 30.94538462  402.290000 29.1900000
#> 45 env1      P-E            T2M_MIN 22.07923077  287.030000 21.6900000
#> 46 env1      P-E                TH1 76.76629085  997.961781 75.3024542
#> 47 env1      P-E                TH2 75.28492308  978.704000 74.0364000
#> 48 env1      P-E                VPD  0.87655065   11.395159  0.6340613
#> 49 env1      P-E               WS2M  0.42461538    5.520000  0.3800000
#> 50 env1      P-E                  n  4.37569231   56.884000  2.0870000
#>      prob_0.5  prob_0.75
#> 1  36.2450000 36.6275000
#> 2   8.9400000 11.4650000
#> 3   8.3150000 17.4975000
#> 4  18.0900000 23.8250000
#> 5   1.2650000  1.5475000
#> 6   0.0400000  0.0500000
#> 7   9.9250000 12.3850000
#> 8   0.0000000  0.0000000
#> 9   0.7650000  0.7775000
#> 10  0.8000000  0.8350000
#> 11 12.5075000 12.5552500
#> 12  0.3540806  0.4347733
#> 13  5.1050000 15.8875000
#> 14 18.5250000 18.8675000
#> 15 86.8000000 88.8650000
#> 16 39.5715000 39.7617500
#> 17 25.7900000 26.5275000
#> 18 23.2550000 23.5325000
#> 19 30.2250000 31.0300000
#> 20 22.3900000 22.7000000
#> 21 76.9598732 77.8659916
#> 22 75.4212000 76.1458000
#> 23  0.6453706  0.7696406
#> 24  0.5350000  0.5775000
#> 25  2.6315000  5.5290000
#> 26 35.9800000 36.2200000
#> 27 10.8000000 11.5900000
#> 28 15.3500000 17.7300000
#> 29 22.3600000 24.1700000
#> 30  1.5100000  1.6100000
#> 31  0.0500000  0.0500000
#> 32  8.0600000 10.5300000
#> 33  0.0000000  0.0000000
#> 34  0.7100000  0.7500000
#> 35  0.7500000  0.8000000
#> 36 12.6730000 12.7010000
#> 37  0.4049855  0.4432030
#> 38  8.1700000 12.1000000
#> 39 17.5100000 18.0400000
#> 40 80.1000000 85.5400000
#> 41 40.1170000 40.1800000
#> 42 26.3800000 26.9800000
#> 43 22.3500000 22.8200000
#> 44 31.2900000 32.8300000
#> 45 22.1800000 22.3200000
#> 46 77.1698742 77.9281685
#> 47 75.6052000 76.2984000
#> 48  0.9479085  1.1868980
#> 49  0.4300000  0.4600000
#> 50  4.8500000  5.5970000
# }
```
