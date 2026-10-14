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
#> start.day is NULL; matched as 2024-01-18
#> end.day is NULL; matched as 2024-02-17
#> ---------------------------------------------------------------
#> get_weather() - daily weather from NASA POWER
#> https://docs.ropensci.org/nasapower
#> ---------------------------------------------------------------
#> Start .............................. Wed Oct 14 03:16:57 2026
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
#> Done ............................... Wed Oct 14 03:16:58 2026
#> Total time ......................... 0.89s
summaryWTH(env.data, env.id = "env", statistic = "mean")
#> ---------------------------------------------------------------
#> summaryWTH -- summarises weather by environment and interval
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: No valid data for: FROST_DAYS. Their summaries are NA/NaN; consider dropping them via 'var.id'.
#>     env           variable        mean
#> 1  env1                T2M 26.31645161
#> 2  env1            T2M_MAX 31.01806452
#> 3  env1            T2M_MIN 22.43258065
#> 4  env1             T2MDEW 22.93387097
#> 5  env1  ALLSKY_SFC_LW_DWN 36.18354839
#> 6  env1  ALLSKY_SFC_SW_DWN 21.17741935
#> 7  env1  ALLSKY_SFC_SW_DNI 13.59032258
#> 8  env1 ALLSKY_SFC_PAR_TOT 10.25645161
#> 9  env1     ALLSKY_SFC_UVA  1.41935484
#> 10 env1     ALLSKY_SFC_UVB  0.04483871
#> 11 env1            PRECTOT  8.74225806
#> 12 env1            EVPTRNS  9.03000000
#> 13 env1               QV2M 18.18451613
#> 14 env1               RH2M 83.30096774
#> 15 env1           GWETROOT  0.72677419
#> 16 env1            GWETTOP  0.77096774
#> 17 env1         FROST_DAYS          NA
#> 18 env1               WS2M  0.45516129
#> 19 env1                VPD  0.81727393
#> 20 env1                  N 12.65741935
#> 21 env1                RTA 40.03732258
#> 22 env1                  n  4.29764516
#> 23 env1                TH1 77.37432518
#> 24 env1                TH2 75.77264516
#> 25 env1           PAR_TEMP  0.38885407

summaryWTH(env.data, env.id = "env", by.interval = TRUE,
           time.window  = c(0, 14, 35, 60, 90, 120),
           names.window = c("P-E", "E-V1", "V1-V4", "V4-VT", "VT-GF", "GF-PM"))
#> ---------------------------------------------------------------
#> summaryWTH -- summarises weather by environment and interval
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: No valid data for: FROST_DAYS. Their summaries are NA/NaN; consider dropping them via 'var.id'.
#>     env interval           variable        mean         sum  prob_0.25
#> 1  env1     E-V1  ALLSKY_SFC_LW_DWN 36.08833333  649.590000 35.9650000
#> 2  env1     E-V1 ALLSKY_SFC_PAR_TOT  9.65500000  173.790000  7.9225000
#> 3  env1     E-V1  ALLSKY_SFC_SW_DNI 11.46833333  206.430000  4.6500000
#> 4  env1     E-V1  ALLSKY_SFC_SW_DWN 19.86166667  357.510000 15.9875000
#> 5  env1     E-V1     ALLSKY_SFC_UVA  1.33944444   24.110000  1.1300000
#> 6  env1     E-V1     ALLSKY_SFC_UVB  0.04333333    0.780000  0.0400000
#> 7  env1     E-V1            EVPTRNS  9.42055556  169.570000  7.6625000
#> 8  env1     E-V1         FROST_DAYS          NA          NA         NA
#> 9  env1     E-V1           GWETROOT  0.74944444   13.490000  0.7300000
#> 10 env1     E-V1            GWETTOP  0.80166667   14.430000  0.7625000
#> 11 env1     E-V1                  N 12.59461111  226.703000 12.5497500
#> 12 env1     E-V1           PAR_TEMP  0.37183143    6.692966  0.3214373
#> 13 env1     E-V1            PRECTOT 10.45666667  188.220000  4.5575000
#> 14 env1     E-V1               QV2M 18.16166667  326.910000 17.8125000
#> 15 env1     E-V1               RH2M 84.98666667 1529.760000 82.6500000
#> 16 env1     E-V1                RTA 39.88094444  717.857000 39.7412500
#> 17 env1     E-V1                T2M 25.90500000  466.290000 24.9650000
#> 18 env1     E-V1             T2MDEW 22.92944444  412.730000 22.6375000
#> 19 env1     E-V1            T2M_MAX 30.42277778  547.610000 29.1075000
#> 20 env1     E-V1            T2M_MIN 22.16388889  398.950000 21.7325000
#> 21 env1     E-V1                TH1 76.88914863 1384.004675 75.7643346
#> 22 env1     E-V1                TH2 75.35960000 1356.472800 74.3871000
#> 23 env1     E-V1                VPD  0.72381251   13.028625  0.6224861
#> 24 env1     E-V1               WS2M  0.48222222    8.680000  0.4525000
#> 25 env1     E-V1                  n  3.62211111   65.198000  1.4682500
#> 26 env1      P-E  ALLSKY_SFC_LW_DWN 36.31538462  472.100000 35.9000000
#> 27 env1      P-E ALLSKY_SFC_PAR_TOT 11.08923077  144.160000  9.8400000
#> 28 env1      P-E  ALLSKY_SFC_SW_DNI 16.52846154  214.870000 10.9200000
#> 29 env1      P-E  ALLSKY_SFC_SW_DWN 22.99923077  298.990000 19.8000000
#> 30 env1      P-E     ALLSKY_SFC_UVA  1.53000000   19.890000  1.3700000
#> 31 env1      P-E     ALLSKY_SFC_UVB  0.04692308    0.610000  0.0400000
#> 32 env1      P-E            EVPTRNS  8.48923077  110.360000  6.8400000
#> 33 env1      P-E         FROST_DAYS          NA          NA         NA
#> 34 env1      P-E           GWETROOT  0.69538462    9.040000  0.6800000
#> 35 env1      P-E            GWETTOP  0.72846154    9.470000  0.7100000
#> 36 env1      P-E                  N 12.74438462  165.677000 12.7200000
#> 37 env1      P-E           PAR_TEMP  0.41242389    5.361511  0.3693694
#> 38 env1      P-E            PRECTOT  6.36846154   82.790000  1.8700000
#> 39 env1      P-E               QV2M 18.21615385  236.810000 17.4600000
#> 40 env1      P-E               RH2M 80.96692308 1052.570000 78.2500000
#> 41 env1      P-E                RTA 40.25384615  523.300000 40.2160000
#> 42 env1      P-E                T2M 26.88615385  349.520000 26.6800000
#> 43 env1      P-E             T2MDEW 22.94000000  298.220000 22.2600000
#> 44 env1      P-E            T2M_MAX 31.84230769  413.950000 31.0300000
#> 45 env1      P-E            T2M_MIN 22.80461538  296.460000 22.2400000
#> 46 env1      P-E                TH1 78.04610810 1014.599405 77.6140656
#> 47 env1      P-E                TH2 76.34455385  992.479200 76.0212000
#> 48 env1      P-E                VPD  0.94668205   12.306867  0.7845527
#> 49 env1      P-E               WS2M  0.41769231    5.430000  0.3500000
#> 50 env1      P-E                  n  5.23300000   68.029000  3.4570000
#>      prob_0.5  prob_0.75
#> 1  36.0300000 36.2425000
#> 2  10.1600000 11.4650000
#> 3  12.9500000 17.3550000
#> 4  20.8950000 23.9075000
#> 5   1.4350000  1.5475000
#> 6   0.0500000  0.0500000
#> 7  10.2200000 12.2875000
#> 8          NA         NA
#> 9   0.7550000  0.7700000
#> 10  0.8000000  0.8275000
#> 11 12.5960000 12.6405000
#> 12  0.3967668  0.4347733
#> 13  9.2350000 13.6050000
#> 14 18.2000000 18.5675000
#> 15 84.7300000 86.5225000
#> 16 39.9030000 40.0350000
#> 17 26.2550000 26.6900000
#> 18 22.9500000 23.3100000
#> 19 30.8550000 31.8725000
#> 20 22.1400000 22.5250000
#> 21 77.0831490 78.0306901
#> 22 75.5802000 76.2768000
#> 23  0.7055348  0.8454940
#> 24  0.4700000  0.5225000
#> 25  4.0880000  5.4795000
#> 26 36.3700000 36.8900000
#> 27 11.4300000 12.5600000
#> 28 16.3600000 22.3900000
#> 29 23.9300000 26.2400000
#> 30  1.5600000  1.7200000
#> 31  0.0500000  0.0500000
#> 32  7.8400000 10.8900000
#> 33         NA         NA
#> 34  0.6900000  0.7100000
#> 35  0.7200000  0.7500000
#> 36 12.7460000 12.7700000
#> 37  0.4209945  0.4590462
#> 38  4.2600000  8.4800000
#> 39 18.1100000 18.9300000
#> 40 80.3700000 84.2100000
#> 41 40.2610000 40.2970000
#> 42 26.7500000 27.0200000
#> 43 22.8900000 23.6200000
#> 44 31.3800000 32.9600000
#> 45 22.8800000 23.3500000
#> 46 78.0612210 78.3335893
#> 47 76.3636000 76.5360000
#> 48  0.9479085  1.1329342
#> 49  0.4200000  0.4500000
#> 50  5.1840000  7.0980000
# }
```
