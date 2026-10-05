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
if (FALSE) { # \dontrun{
wb <- water_balance(wth, soil)

## Season-level, one row per environment
summary_water_balance(wb)

## By maize growth stage
summary_water_balance(wb, by.interval = TRUE,
                      time.window  = c(0, 14, 35, 60, 90, 120),
                      names.window = c("P-E", "E-V1", "V1-V4",
                                       "V4-VT", "VT-GF", "GF-PM"))

## Default 10-day windows when no breakpoints are given
summary_water_balance(wb, by.interval = TRUE)

## Rank environments by flowering-period stress
s <- summary_water_balance(wb, by.interval = TRUE,
                           time.window  = c(0, 60, 90, 120),
                           names.window = c("veg", "flower", "fill"))
s[s$interval == "flower", ][order(s$mean_Ks[s$interval == "flower"]), ]
} # }
```
