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
if (FALSE) { # \dontrun{
data("maizeWTH")

## 1. Default maize staging on cumulative GDD
ph <- env_phenology(maizeWTH, crop = "maize")
table(ph$env, ph$stage)

## 2. Photothermal time, daylength taken from param_radiation()'s N column
ph2 <- env_phenology(processWTH(maizeWTH), crop = "maize",
                     photoperiod = TRUE, daylength.id = "N")

## 3. Custom stages and the Baskerville-Emin method
ph3 <- env_phenology(maizeWTH,
                     stages = c(veg = 500, flow = 1200, fill = 2200),
                     Tbase = 8, Tupper = 32, method = "baskerville")

## 4. Hand the stage factor to the existing covariable machinery.
##    This is the point of the function: W columns now align on
##    development instead of on calendar date.
W <- W_matrix(env.data = ph, var.id = c("T2M", "PRECTOT"),
              by.interval = TRUE, time.window = NULL)
} # }
```
