# Thermal-Time Phenology: Map Daily Weather onto Developmental Stages

Accumulates daily thermal time from planting and cuts each environment's
season at crop-specific developmental thresholds, so that every
downstream interval-based routine compares LIKE DEVELOPMENTAL STATES
across sites rather than like calendar dates.

The returned `stage` column is a drop-in replacement for the `interval`
produced by `time.window` in
[`W_matrix()`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md),
[`summaryWTH()`](https://gcostaneto.github.io/EnvRtype/reference/summaryWTH.md)
and
[`env_typing()`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md).

## Usage

``` r
env_phenology(
  env.data,
  env.id = "env",
  day.id = "daysFromStart",
  Tmax = "T2M_MAX",
  Tmin = "T2M_MIN",
  crop = NULL,
  stages = NULL,
  Tbase = NULL,
  Tupper = NULL,
  method = c("gdd", "baskerville", "trapezoid"),
  Tmaxdev = 45,
  photoperiod = FALSE,
  lat = "LAT",
  daylength.id = NULL,
  p.opt = NULL,
  p.sens = 0.02,
  day.type = NULL,
  planting.id = NULL,
  merge = TRUE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. Daily weather, typically a
  [`get_weather()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md)
  or
  [`processWTH()`](https://gcostaneto.github.io/EnvRtype/reference/processWTH.md)
  output.

- env.id:

  character. Environment-id column. Default `"env"`.

- day.id:

  character. Day-counter column. Default `"daysFromStart"`.

- Tmax, Tmin:

  character. Daily maximum/minimum temperature columns. Defaults
  `"T2M_MAX"` / `"T2M_MIN"`.

- crop:

  character. Name of a built-in template, e.g. `"maize"`, `"wheat"`,
  `"soybean"`, `"rice"`, `"cotton"`, `"canola"`. Supplies default
  `stages`, `Tbase`, `Tupper`, `p.opt` and `day.type`. Any explicitly
  supplied argument wins. Call
  [`show_phenology()`](https://gcostaneto.github.io/EnvRtype/reference/phenology_templates.md)
  for the full list of the 20 templates and their provenance.

- stages:

  named numeric. Cumulative thermal time (C d) at the END of each stage,
  strictly increasing, e.g.
  `c(VE = 125, V6 = 475, VT = 1135, R6 = 2700)`. Overrides `crop`.

- Tbase, Tupper:

  numeric. Cardinal temperatures for thermal time. Override the crop
  template when supplied.

- method:

  character. Thermal-time method: `"gdd"` (default, McMaster & Wilhelm
  with upper cutoff), `"baskerville"` (Baskerville-Emin sine
  integration), or `"trapezoid"` (linear supra-optimal decline to
  `Tmaxdev`).

- Tmaxdev:

  numeric. Absolute developmental maximum, used by
  `method = "trapezoid"`. Default 45.

- photoperiod:

  boolean. If `TRUE`, multiply daily thermal time by a photoperiod
  factor, producing photothermal time (PTT).

- lat:

  character or numeric. Latitude column name or a numeric vector, used
  to derive daylength when `photoperiod = TRUE` and no daylength column
  exists. Default `"LAT"`.

- daylength.id:

  character. Existing daylength column (h), e.g. `"N"` from
  [`param_radiation()`](https://gcostaneto.github.io/EnvRtype/reference/param_radiation.md).
  Used in preference to deriving it.

- p.opt, p.sens:

  numeric. Optimum daylength (h) and sensitivity (fractional slowdown
  per hour of deviation). When `p.opt` is `NULL` (default) the crop
  template's value is used – 12.5 h for short-day maize, 16 h for
  long-day wheat, and so on.

- day.type:

  character. `"short"`, `"long"` or `"neutral"`. `NULL` (default) takes
  the crop template's value. A `"neutral"` crop silently keeps a
  photoperiod factor of 1 even when `photoperiod = TRUE`.

- planting.id:

  character. Optional column giving planting date, used to reset
  accumulation. If `NULL`, accumulation starts at the first record of
  each environment.

- merge:

  boolean. If `TRUE` (default) return `env.data` with the new columns
  appended; if `FALSE` return only the new columns.

- verbose:

  boolean. Print a progress banner.

## Value

A data.frame with these columns added:

- dTT:

  daily thermal time (C d), or photothermal time if `photoperiod = TRUE`

- cumTT:

  cumulative thermal time from planting

- stage:

  ordered factor, the developmental stage of that day

- dstage:

  relative position within the stage, in \[0, 1)

- stage_day:

  day counter within the stage (1-based)

The attribute `"phenology"` carries the stage table, the method and the
cardinal temperatures. Environments whose season ends before the last
threshold get `NA` for the stages never reached, and a warning names
them – silently recycling the final stage would fabricate development
that the weather does not support.

## Details

**Which method.** `"gdd"` is the convention and what
[`param_temperature()`](https://gcostaneto.github.io/EnvRtype/reference/param_temperature.md)
already computes; use it for comparability. `"baskerville"` is
materially better where daily range straddles `Tbase` – high-elevation
or early-spring sites – because the simple mean then understates
accumulation. `"trapezoid"` is the only option that penalises
supra-optimal heat, relevant for heat-stress work.

**Photothermal time.** Maize, sorghum and rice are quantitative
short-day plants: long days delay flowering. With `photoperiod = TRUE`
the daily increment is scaled by a linear daylength factor, which keeps
tropical and temperate plantings on a common developmental scale. The
default `p.sens = 0.02` is deliberately mild; calibrate it before
relying on the absolute stage dates.

**Caveat.** These thresholds are cultivar-generic. Maturity group alone
moves maize R6 by 300-500 C d. Treat the stage boundaries as a
consistent relative scale, not as a validated phenology prediction,
unless you have calibrated `stages` against observed dates.

## References

McMaster, G.S. & Wilhelm, W.W. (1997) Growing degree-days: one equation,
two interpretations. *Agric. For. Meteorol.* 87(4), 291-300.

Baskerville, G.L. & Emin, P. (1969) Rapid estimation of heat
accumulation from maximum and minimum temperatures. *Ecology* 50(3),
514-517.

Abendroth, L.J. et al. (2011) *Corn Growth and Development*. PMR 1009,
Iowa State University Extension.

## See also

[`param_temperature`](https://gcostaneto.github.io/EnvRtype/reference/param_temperature.md),
[`W_matrix`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md),
[`env_risk_profile`](https://gcostaneto.github.io/EnvRtype/reference/env_risk_profile.md),
[`summaryWTH`](https://gcostaneto.github.io/EnvRtype/reference/summaryWTH.md)

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
data("maizeWTH")

## 1. Default maize staging on cumulative GDD
ph <- env_phenology(maizeWTH, crop = "maize")
#> ---------------------------------------------------------------
#> env_phenology -- maps weather onto thermal-time stages
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing daily thermal time (method = gdd)
#>   - assigning developmental stages
#> Warning: 5 environment(s) did not accumulate enough thermal time to reach 'R6': IP, NM, PM, SE, SO. Days beyond the last reached stage are NA. Extend the series or lower the stage thresholds.
#> ---------------------------------------------------------------
#> Thermal-time stages (GDD, Tbase = 10 C, Tupper = 30 C)
#>   VE      <=     125 C d       40 day(s)
#>   V6      <=     475 C d      123 day(s)
#>   V12     <=     740 C d       94 day(s)
#>   VT      <=    1135 C d      152 day(s)
#>   R1      <=    1250 C d       47 day(s)
#>   R3      <=    1660 C d      161 day(s)
#>   R5      <=    2110 C d      107 day(s)
#>   R6      <=    2700 C d       32 day(s)
#> ---------------------------------------------------------------
table(ph$env, ph$stage)
#>     
#>      VE V6 V12 VT R1 R3 R5 R6
#>   NM  8 23  18 26  8 28 25  0
#>   SO  7 22  16 25  7 26 29 22
#>   PM  8 27  20 34 11 36 17  0
#>   IP  8 24  18 27  9 32 36 10
#>   SE  9 27  22 40 12 39  0  0

## 2. Photothermal time, daylength taken from param_radiation()'s N column
ph2 <- env_phenology(processWTH(maizeWTH), crop = "maize",
                     photoperiod = TRUE, daylength.id = "N")
#> ---------------------------------------------------------------
#> env_phenology -- maps weather onto thermal-time stages
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> ---------------------------------------------------------------
#> processWTH -- derives all agro-meteorological parameters
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Dropping existing column(s) before recomputation: n, N, RTA, VPD, SPV, ETP, PETP, GDD, FRUE, T2M_RANGE
#>   - radiation parameters computed
#>   - atmospheric-demand parameters computed
#>   - temperature-response parameters computed
#>   - computing daily thermal time (method = gdd)
#>   - thermal time modulated by photoperiod (PTT)
#>   - assigning developmental stages
#> Warning: 5 environment(s) did not accumulate enough thermal time to reach 'R6': IP, NM, PM, SE, SO. Days beyond the last reached stage are NA. Extend the series or lower the stage thresholds.
#> ---------------------------------------------------------------
#> Thermal-time stages (PTT, Tbase = 10 C, Tupper = 30 C)
#>   VE      <=     125 C d       40 day(s)
#>   V6      <=     475 C d      123 day(s)
#>   V12     <=     740 C d       94 day(s)
#>   VT      <=    1135 C d      152 day(s)
#>   R1      <=    1250 C d       47 day(s)
#>   R3      <=    1660 C d      161 day(s)
#>   R5      <=    2110 C d      107 day(s)
#>   R6      <=    2700 C d       32 day(s)
#> ---------------------------------------------------------------

## 3. Custom stages and the Baskerville-Emin method
ph3 <- env_phenology(maizeWTH,
                     stages = c(veg = 500, flow = 1200, fill = 2200),
                     Tbase = 8, Tupper = 32, method = "baskerville")
#> ---------------------------------------------------------------
#> env_phenology -- maps weather onto thermal-time stages
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing daily thermal time (method = baskerville)
#>   - assigning developmental stages
#> Warning: 2 environment(s) did not accumulate enough thermal time to reach 'fill': PM, SE. Days beyond the last reached stage are NA. Extend the series or lower the stage thresholds.
#> ---------------------------------------------------------------
#> Thermal-time stages (GDD, Tbase = 8 C, Tupper = 32 C)
#>   veg     <=     500 C d      149 day(s)
#>   flow    <=    1200 C d      223 day(s)
#>   fill    <=    2200 C d      314 day(s)
#> ---------------------------------------------------------------

## 4. Hand the stage factor to the existing covariable machinery.
##    This is the point of the function: W columns now align on
##    development instead of on calendar date.
W <- W_matrix(env.data = ph, var.id = c("T2M", "PRECTOT"),
              by.interval = TRUE, time.window = NULL)
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W
# }
```
