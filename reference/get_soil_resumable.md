# Crash-Resilient Wrapper for get_soil()

Downloads SoilGrids properties one environment at a time with the same
checkpointing scheme as `get_weather_resumable`: each result is written
to disk immediately and recorded in an append-only log, so a crashed run
resumes with `restart_from_log`.

SoilGrids is rate-limited and slower than NASA POWER, which makes long
soil jobs the most frequent victims of a mid-run failure.

## Usage

``` r
get_soil_resumable(
  env.id = NULL,
  lat = NULL,
  lon = NULL,
  variables.names = "clay",
  depths = c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"),
  stat = c("mean", "q_5", "q_50", "q_95", "uncertainty"),
  dir.path = NULL,
  run.id = "soil_run",
  sleep = 20,
  tries = 3L,
  resume = TRUE,
  verbose = TRUE
)
```

## Arguments

- env.id, lat, lon, variables.names, depths, stat:

  Passed to `get_soil`.

- dir.path:

  character. Parent directory. Default
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- run.id:

  character. Run directory name. Default `"soil_run"`.

- sleep:

  numeric. Seconds between environments. Default 20, matching
  `get_soil`; lowering it much invites throttling by ISRIC.

- tries:

  integer. Retry attempts per environment. Default 3.

- resume:

  logical. Skip environments already logged `ok`. Default `TRUE`.

- verbose:

  logical. Print progress. Default `TRUE`.

## Value

A `data.frame` of soil properties (wide by default), with attributes
`run.dir`, `n_ok`, `n_failed` and `failed`.

## Details

Results are stored **long** on disk, one row per property-depth, and
reshaped to wide only at assembly. Storing long is what makes resumption
safe: parts fetched at different times still stack correctly, whereas
pre-widened parts could have mismatched columns.

See `get_weather_resumable` for the full checkpoint design.

## See also

`restart_from_log`, `read_progress_log`, `get_soil`

## Examples

``` r
if (FALSE) { # \dontrun{
grid <- data.frame(env = paste0("P", 1:200),
                   lat = runif(200, -13, -12),
                   lon = runif(200, -56, -55))

soil <- get_soil_resumable(
  env.id = grid$env, lat = grid$lat, lon = grid$lon,
  variables.names = c("clay", "sand", "silt", "soc", "phh2o"),
  run.id = "zoning_grid")

## After a crash, in a new session:
soil <- restart_from_log(run.id = "zoning_grid")

## Points in the ocean never resolve -- keep the successes
soil <- restart_from_log(run.id = "zoning_grid", retry.failed = FALSE)
} # }
```
