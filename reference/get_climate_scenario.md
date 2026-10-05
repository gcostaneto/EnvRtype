# Future-Climate Envirotyping from CMIP6 Downscaled Scenarios

Retrieves downscaled CMIP6 monthly climate (WorldClim 2.1) for a set of
sites and converts it into a DAILY weather series compatible with
[`processWTH()`](https://gcostaneto.github.io/EnvRtype/reference/processWTH.md),
[`env_phenology()`](https://gcostaneto.github.io/EnvRtype/reference/env_phenology.md)
and
[`water_balance()`](https://gcostaneto.github.io/EnvRtype/reference/water_balance.md),
by applying the monthly change signal to an observed baseline series
(delta-change downscaling).

## Usage

``` r
get_climate_scenario(
  env.id = NULL,
  lat = NULL,
  lon = NULL,
  baseline = NULL,
  model = "MPI-ESM1-2-HR",
  models = c("IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"),
  scenario = "ssp245",
  period = "2041-2060",
  resolution = c("2.5m", "5m", "10m", "30s"),
  method = c("delta", "raw"),
  dir.path = NULL,
  baseline.id = "env",
  baseline.date = "YYYYMMDD",
  precip.scale = c("multiplicative", "additive"),
  verbose = TRUE
)
```

## Arguments

- env.id:

  vector (character). Site identifiers.

- lat, lon:

  vector (numeric). Site coordinates.

- baseline:

  data.frame. Observed daily weather for the same sites, from
  [`get_weather()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md).
  Supplies the day-to-day variability that monthly projections cannot.
  REQUIRED for `method = "delta"`.

- model:

  character. CMIP6 GCM name, e.g. `"MPI-ESM1-2-HR"`. One of `.CS_GCMS`.
  Use `"ensemble"` to average several.

- models:

  vector (character). Used when `model = "ensemble"`. Default is a
  four-model spread across climate sensitivities.

- scenario:

  character. `"ssp126"`, `"ssp245"`, `"ssp370"` or `"ssp585"`.

- period:

  character. `"2021-2040"`, `"2041-2060"`, `"2061-2080"` or
  `"2081-2100"`.

- resolution:

  character. WorldClim tile resolution: `"2.5m"` (default), `"5m"`,
  `"10m"` or `"30s"`. Finer is a much larger download. Resolution does
  NOT affect which model x scenario combinations exist: archive gaps are
  per GCM x SSP and are identical at every resolution and period. See
  `.CS_GAPS`.

- method:

  character. `"delta"` (default) perturbs the baseline; `"raw"` returns
  the monthly projected values without downscaling.

- dir.path:

  character. Cache directory for downloaded rasters. Defaults to
  `tools::R_user_dir("et", "cache")`.

- baseline.id, baseline.date:

  character. Columns in `baseline`.

- precip.scale:

  character. `"multiplicative"` (default) or `"additive"` application of
  the precipitation delta.

- verbose:

  boolean. Print progress.

## Value

For `method = "delta"`, a data.frame in
[`get_weather()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md)
layout (`env`, `YYYYMMDD`, `daysFromStart`, `T2M`, `T2M_MAX`, `T2M_MIN`,
`PRECTOT`, ...) with the scenario signal applied, plus `scenario`,
`model` and `period` columns. For `method = "raw"`, a monthly data.frame
of projected `tmin`, `tmax` and `prec` per site. The attribute
`"deltas"` holds the monthly change factors that were applied, so the
perturbation is auditable.

## Details

**What delta-change does and does not do.** The monthly mean change from
the GCM is added to (temperature) or multiplied onto (precipitation) the
observed daily series. The result therefore keeps the OBSERVED
day-to-day sequence, wet-day frequency and spell structure, and shifts
only the monthly means. This is the standard, defensible way to drive a
crop model from monthly projections (Hawkins et al. 2013), but note the
consequence: **changes in variability, in dry-spell length, and in
extreme-event frequency are NOT represented**. A scenario that in
reality brings the same total rain in fewer, heavier events will look
identical to one that brings it evenly. Heat-stress day counts under
delta-change are driven purely by the mean shift, so they are a lower
bound.

If extremes matter to your question, you need daily GCM output and a
quantile-mapping bias correction, not this function.

**Model choice is not innocuous.** Across CMIP6, equilibrium climate
sensitivity ranges roughly 1.8-5.6 K, and several high-sensitivity
models are now regarded as implausibly hot. A single-model result is one
draw from a wide distribution. Prefer `model = "ensemble"` and report
the spread; a projection whose sign flips across models is not a
finding.

**Network access.** WorldClim is an external host. In restricted
environments the download will fail; pre-populate `dir.path` with the
GeoTIFFs and the function will use the cache.

## References

Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2. *Int. J. Climatol.*
37(12), 4302-4315.

Hawkins, E. et al. (2013) Calibration and bias correction of climate
projections for crop modelling. *Agric. For. Meteorol.* 170, 19-31.

O'Neill, B.C. et al. (2016) The Scenario Model Intercomparison Project
(ScenarioMIP) for CMIP6. *Geosci. Model Dev.* 9, 3461-3482.

## See also

[`get_weather`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md),
[`project_risk`](https://gcostaneto.github.io/EnvRtype/reference/project_risk.md),
[`env_phenology`](https://gcostaneto.github.io/EnvRtype/reference/env_phenology.md)

## Author

Germano Costa Neto

## Examples

``` r
if (FALSE) { # \dontrun{
sites <- data.frame(env = c("PIRA", "SETE"),
                    lat = c(-22.7, -19.4), lon = c(-47.6, -44.2))

base <- get_weather(env.id = sites$env, lat = sites$lat, lon = sites$lon,
                    start.day = "2015-10-01", end.day = "2016-03-31")

## Mid-century, intermediate scenario, 4-model ensemble.
## The default 'models' set covers all four SSPs.
fut <- get_climate_scenario(env.id = sites$env, lat = sites$lat,
                            lon = sites$lon, baseline = base,
                            model = "ensemble", scenario = "ssp245",
                            period = "2041-2060")

## Which GCMs fully cover a scenario? (GFDL-ESM4 has no ssp245.)
.cs_models_for("ssp245")

## GFDL-ESM4 is fine for ssp126 or ssp370, just not ssp245/ssp585:
fut126 <- get_climate_scenario(env.id = sites$env, lat = sites$lat,
                               lon = sites$lon, baseline = base,
                               model = "GFDL-ESM4", scenario = "ssp126",
                               period = "2041-2060")

## Drop straight into the existing pipeline
fut <- processWTH(fut)
ph  <- env_phenology(fut, crop = "maize")

## Compare thermal time now vs mid-century
tapply(ph$dTT, ph$env, sum)
} # }
```
