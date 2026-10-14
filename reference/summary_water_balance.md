# Summarise Water-Balance Output by Environment or Growth Stage

Condenses the daily output of `water_balance` into one row per
environment, or per environment and growth stage, producing the
season-level indices most often used as environmental covariates in
genotype-by-environment models.

## Usage

``` r
summary_water_balance(
  wb,
  env.id = "env",
  days.id = "daysFromStart",
  by.interval = FALSE,
  time.window = NULL,
  names.window = NULL,
  verbose = TRUE
)
```

## Arguments

- wb:

  data.frame. Output of `water_balance`.

- env.id:

  character. Environment-id column. Default `"env"`.

- days.id:

  character. Day counter column. Default `"daysFromStart"`.

- by.interval:

  logical. Summarise within growth stages as well as by environment.
  Default `FALSE`.

- time.window:

  numeric. Stage breakpoints in days after planting. When `NULL`, 10-day
  windows are used.

- names.window:

  character. Stage names, one per interval. Ignored if the length does
  not match the number of intervals.

- verbose:

  boolean. If `TRUE` (default) prints a progress banner.

## Value

A data.frame with, per group: `n_days`, `mean_Ks`, `min_Ks`,
`total_deficit`, `n_stress_days`, `max_consecutive_stress`, `total_ETc`,
`total_ETa`, `ETa_ETc_ratio`, `total_rain`, `total_irrigation`,
`total_drainage`, `mean_ASW`, `mean_depletion`.

## Details

`ETa_ETc_ratio` is the relative transpiration deficit, the quantity
entering most water-limited yield models – notably the FAO Ky framework
of Doorenbos & Kassam (1979), where relative yield loss is proportional
to \\1 - ET_a/ET_c\\.

`max_consecutive_stress` distinguishes one damaging spell from the same
number of scattered days, which a mean cannot. Ten consecutive stressed
days at flowering is a different event from ten isolated days across the
season, and yield responds very differently to the two.

**Stage timing matters more than season totals.** Water deficit during
flowering and grain fill costs disproportionately more yield than the
same deficit during vegetative growth (Steduto et al. 2012). Prefer
`by.interval = TRUE` with breakpoints matched to your crop's phenology
over a single season-long mean.

## References

Doorenbos, J. & Kassam, A.H. (1979) Yield response to water. *FAO
Irrigation and Drainage Paper 33*. FAO, Rome.

Steduto, P., Hsiao, T.C., Fereres, E. & Raes, D. (2012) Crop yield
response to water. *FAO Irrigation and Drainage Paper 66*. FAO, Rome.

Costa-Neto, G., Fritsche-Neto, R. & Crossa, J. (2021) Nonlinear kernels,
dominance, and envirotyping data increase the accuracy of genome-based
prediction in multi-environment trials. *Heredity* 126, 92-106.
[doi:10.1038/s41437-020-00353-1](https://doi.org/10.1038/s41437-020-00353-1)

## See also

`water_balance`, `summaryWTH`

## Examples

``` r
# \donttest{
## Synthetic daily weather and a uniform soil profile (mapped units).
wth <- expand.grid(daysFromStart = 1:120, env = c("E1", "E2"))
wth$PRECTOT <- ifelse(wth$daysFromStart %% 10 == 0, 25, 0)
wth$ETP <- 5
soil <- data.frame(env = c("E1", "E2"))
for (layer in c("0_5cm", "5_15cm", "15_30cm", "30_60cm", "60_100cm")) {
  soil[[paste0("wv0033|", layer)]] <- 300
  soil[[paste0("wv1500|", layer)]] <- 150
}
wb <- water_balance(wth, soil)
#> ---------------------------------------------------------------
#> water_balance -- runs a FAO-56 daily soil water balance
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Using 'ETP' as reference evapotranspiration.
#> ---------------------------------------------------------------
#> water_balance() - FAO-56 soil water balance
#> Allen et al. (1998), FAO Irrigation and Drainage Paper 56
#> ---------------------------------------------------------------
#> Environments ....................... 2
#> Rooting depth ...................... 300 -> 1000 mm by 60 DAP
#> Depletion fraction p ............... 0.55 (adjusted daily for ETc)
#>   [1/2] E1 | TAW 47-150 mm | stressed days 79/120
#>   [2/2] E2 | TAW 47-150 mm | stressed days 79/120
#> Stress days: none=82  mild=18  moderate=90  severe=50
#> Done.

## Season-level, one row per environment
summary_water_balance(wb)
#> ---------------------------------------------------------------
#> summary_water_balance -- summarises a water balance by environment/stage
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   env n_days mean_Ks min_Ks total_deficit n_stress_days max_consecutive_stress
#> 1  E1    120  0.7295 0.3119       167.684            79                     77
#> 2  E2    120  0.7295 0.3119       167.684            79                     77
#>   total_ETc total_ETa ETa_ETc_ratio total_rain total_irrigation total_drainage
#> 1   565.875   398.191        0.7037        300                0          0.855
#> 2   565.875   398.191        0.7037        300                0          0.855
#>   mean_ASW mean_depletion
#> 1  44.4743         0.5815
#> 2  44.4743         0.5815

## By maize growth stage
summary_water_balance(wb, by.interval = TRUE,
                      time.window  = c(0, 14, 35, 60, 90, 120),
                      names.window = c("P-E", "E-V1", "V1-V4",
                                       "V4-VT", "VT-GF", "GF-PM"))
#> ---------------------------------------------------------------
#> summary_water_balance -- summarises a water balance by environment/stage
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>    env interval n_days mean_Ks min_Ks total_deficit n_stress_days
#> 1   E1      P-E     13  1.0000 1.0000         0.000             0
#> 2   E1     E-V1     21  1.0000 1.0000         0.000             0
#> 3   E1    V1-V4     25  0.8380 0.5785        23.192            18
#> 4   E1    V4-VT     30  0.5341 0.3390        80.374            30
#> 5   E1    VT-GF     30  0.5359 0.3119        62.601            30
#> 6   E1    GF-PM      1  0.4944 0.4944         1.517             1
#> 7   E2      P-E     13  1.0000 1.0000         0.000             0
#> 8   E2     E-V1     21  1.0000 1.0000         0.000             0
#> 9   E2    V1-V4     25  0.8380 0.5785        23.192            18
#> 10  E2    V4-VT     30  0.5341 0.3390        80.374            30
#> 11  E2    VT-GF     30  0.5359 0.3119        62.601            30
#> 12  E2    GF-PM      1  0.4944 0.4944         1.517             1
#>    max_consecutive_stress total_ETc total_ETa ETa_ETc_ratio total_rain
#> 1                       0    33.184    33.184        1.0000         25
#> 2                       0    85.216    85.216        1.0000         50
#> 3                      16   139.350   116.158        0.8336         50
#> 4                      30   172.500    92.126        0.5341         75
#> 5                      30   132.625    70.024        0.5280         75
#> 6                       1     3.000     1.483        0.4943         25
#> 7                       0    33.184    33.184        1.0000         25
#> 8                       0    85.216    85.216        1.0000         50
#> 9                      16   139.350   116.158        0.8336         50
#> 10                     30   172.500    92.126        0.5341         75
#> 11                     30   132.625    70.024        0.5280         75
#> 12                      1     3.000     1.483        0.4943         25
#>    total_irrigation total_drainage mean_ASW mean_depletion
#> 1                 0          0.855  48.1002         0.1595
#> 2                 0          0.000  60.1676         0.3016
#> 3                 0          0.000  49.3071         0.6056
#> 4                 0          0.000  37.8824         0.7475
#> 5                 0          0.000  34.2662         0.7716
#> 6                 0          0.000  50.9550         0.6603
#> 7                 0          0.855  48.1002         0.1595
#> 8                 0          0.000  60.1676         0.3016
#> 9                 0          0.000  49.3071         0.6056
#> 10                0          0.000  37.8824         0.7475
#> 11                0          0.000  34.2662         0.7716
#> 12                0          0.000  50.9550         0.6603

## Default 10-day windows when no breakpoints are given
summary_water_balance(wb, by.interval = TRUE)
#> ---------------------------------------------------------------
#> summary_water_balance -- summarises a water balance by environment/stage
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>    env    interval n_days mean_Ks min_Ks total_deficit n_stress_days
#> 1   E1  Interval_1      9  1.0000 1.0000         0.000             0
#> 2   E1  Interval_2     10  1.0000 1.0000         0.000             0
#> 3   E1  Interval_3     10  1.0000 1.0000         0.000             0
#> 4   E1  Interval_4     10  0.9934 0.9335         0.351             1
#> 5   E1  Interval_5     10  0.8657 0.6605         7.662             7
#> 6   E1  Interval_6     10  0.7360 0.5785        15.179            10
#> 7   E1  Interval_7     10  0.6182 0.4384        21.956            10
#> 8   E1  Interval_8     10  0.5146 0.3691        27.913            10
#> 9   E1  Interval_9     10  0.4695 0.3390        30.505            10
#> 10  E1 Interval_10     10  0.4684 0.3119        28.252            10
#> 11  E1 Interval_11     10  0.5229 0.3359        20.988            10
#> 12  E1 Interval_12     10  0.6163 0.3954        13.361            10
#> 13  E1 Interval_13      1  0.4944 0.4944         1.517             1
#> 14  E2  Interval_1      9  1.0000 1.0000         0.000             0
#> 15  E2  Interval_2     10  1.0000 1.0000         0.000             0
#> 16  E2  Interval_3     10  1.0000 1.0000         0.000             0
#> 17  E2  Interval_4     10  0.9934 0.9335         0.351             1
#> 18  E2  Interval_5     10  0.8657 0.6605         7.662             7
#> 19  E2  Interval_6     10  0.7360 0.5785        15.179            10
#> 20  E2  Interval_7     10  0.6182 0.4384        21.956            10
#> 21  E2  Interval_8     10  0.5146 0.3691        27.913            10
#> 22  E2  Interval_9     10  0.4695 0.3390        30.505            10
#> 23  E2 Interval_10     10  0.4684 0.3119        28.252            10
#> 24  E2 Interval_11     10  0.5229 0.3359        20.988            10
#> 25  E2 Interval_12     10  0.6163 0.3954        13.361            10
#> 26  E2 Interval_13      1  0.4944 0.4944         1.517             1
#>    max_consecutive_stress total_ETc total_ETa ETa_ETc_ratio total_rain
#> 1                       0    21.316    21.316        1.0000          0
#> 2                       0    32.434    32.434        1.0000         25
#> 3                       0    41.100    41.100        1.0000         25
#> 4                       1    49.100    48.749        0.9929         25
#> 5                       6    56.300    48.638        0.8639         25
#> 6                      10    57.500    42.321        0.7360         25
#> 7                      10    57.500    35.544        0.6182         25
#> 8                      10    57.500    29.587        0.5146         25
#> 9                      10    57.500    26.995        0.4695         25
#> 10                     10    53.375    25.122        0.4707         25
#> 11                     10    44.208    23.220        0.5252         25
#> 12                     10    35.042    21.682        0.6187         25
#> 13                      1     3.000     1.483        0.4943         25
#> 14                      0    21.316    21.316        1.0000          0
#> 15                      0    32.434    32.434        1.0000         25
#> 16                      0    41.100    41.100        1.0000         25
#> 17                      1    49.100    48.749        0.9929         25
#> 18                      6    56.300    48.638        0.8639         25
#> 19                     10    57.500    42.321        0.7360         25
#> 20                     10    57.500    35.544        0.6182         25
#> 21                     10    57.500    29.587        0.5146         25
#> 22                     10    57.500    26.995        0.4695         25
#> 23                     10    53.375    25.122        0.4707         25
#> 24                     10    44.208    23.220        0.5252         25
#> 25                     10    35.042    21.682        0.6187         25
#> 26                      1     3.000     1.483        0.4943         25
#>    total_irrigation total_drainage mean_ASW mean_depletion
#> 1                 0          0.000  42.5219         0.2006
#> 2                 0          0.855  56.1251         0.1939
#> 3                 0          0.000  61.3250         0.2953
#> 4                 0          0.000  58.3600         0.4407
#> 5                 0          0.000  50.4072         0.5860
#> 6                 0          0.000  47.6630         0.6582
#> 7                 0          0.000  43.4536         0.7103
#> 8                 0          0.000  36.5898         0.7561
#> 9                 0          0.000  33.6038         0.7760
#> 10                0          0.000  32.6264         0.7825
#> 11                0          0.000  33.7186         0.7752
#> 12                0          0.000  36.4535         0.7570
#> 13                0          0.000  50.9550         0.6603
#> 14                0          0.000  42.5219         0.2006
#> 15                0          0.855  56.1251         0.1939
#> 16                0          0.000  61.3250         0.2953
#> 17                0          0.000  58.3600         0.4407
#> 18                0          0.000  50.4072         0.5860
#> 19                0          0.000  47.6630         0.6582
#> 20                0          0.000  43.4536         0.7103
#> 21                0          0.000  36.5898         0.7561
#> 22                0          0.000  33.6038         0.7760
#> 23                0          0.000  32.6264         0.7825
#> 24                0          0.000  33.7186         0.7752
#> 25                0          0.000  36.4535         0.7570
#> 26                0          0.000  50.9550         0.6603

## Rank environments by flowering-period stress
s <- summary_water_balance(wb, by.interval = TRUE,
                           time.window  = c(0, 60, 90, 120),
                           names.window = c("veg", "flower", "fill"))
#> ---------------------------------------------------------------
#> summary_water_balance -- summarises a water balance by environment/stage
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
s[s$interval == "flower", ][order(s$mean_Ks[s$interval == "flower"]), ]
#>  [1] env                    interval               n_days                
#>  [4] mean_Ks                min_Ks                 total_deficit         
#>  [7] n_stress_days          max_consecutive_stress total_ETc             
#> [10] total_ETa              ETa_ETc_ratio          total_rain            
#> [13] total_irrigation       total_drainage         mean_ASW              
#> [16] mean_depletion        
#> <0 rows> (or 0-length row.names)
# }
```
