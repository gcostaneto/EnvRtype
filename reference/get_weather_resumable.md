# Crash-Resilient Wrapper for get_weather()

Downloads daily weather one environment at a time, writing each result
to disk and appending a line to a progress log as soon as it succeeds.
If the session dies part-way through, `restart_from_log` resumes from
the last completed environment instead of restarting the whole job.

Use this instead of `get_weather` whenever losing the run would hurt –
in practice, more than a few dozen environments.

## Usage

``` r
get_weather_resumable(
  env.id = NULL,
  lat = NULL,
  lon = NULL,
  start.day = NULL,
  end.day = NULL,
  variables.names = NULL,
  dir.path = NULL,
  run.id = "weather_run",
  sleep = 1,
  tries = 3L,
  resume = TRUE,
  verbose = TRUE
)
```

## Arguments

- env.id, lat, lon, start.day, end.day, variables.names:

  Passed to `get_weather`. `lat` and `lon` are required.

- dir.path:

  character. Parent directory for run directories. Default
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- run.id:

  character. Name of this run's directory. Default `"weather_run"`. Use
  a distinct name per job.

- sleep:

  numeric. Seconds between environments. Default 1.

- tries:

  integer. Retry attempts per environment before it is logged as failed.
  Default 3.

- resume:

  logical. If `TRUE` (default), environments already marked `ok` in an
  existing log are skipped.

- verbose:

  logical. Print progress. Default `TRUE`.

## Value

A `data.frame` of daily weather for every environment that succeeded,
with attributes `run.dir`, `n_ok`, `n_failed` and `failed`.

## Details

**Why per-environment files.** Results are written as one `.rds` per
environment under `parts/`, never accumulated only in memory. A process
killed at any point loses at most the environment in flight.

**Why an append-only log.** Each line is appended and flushed
individually. Rewriting one state file every iteration would risk the
whole file to a kill mid-write; appending risks only the last line, and
the reader tolerates a truncated final record.

**Why a fingerprint.** The manifest stores a hash of coordinates, dates
and variables. Resuming with different arguments is refused, because
silently concatenating two different queries produces a table that looks
fine and means nothing.

**Failures do not abort the run.** A failed environment is logged and
the loop continues. Inspect with `read_progress_log`, then retry with
`restart_from_log`.

**Serial by design.** Unlike `get_weather`, this wrapper does not
parallelise. Interleaved writes from several workers make progress
ordering ambiguous, and NASA POWER rate-limits aggressive clients
anyway. The cost is wall-clock time; the benefit is a log you can trust.
For speed on a job you can afford to lose, use
`get_weather(parallel = TRUE)`.

## See also

`restart_from_log`, `read_progress_log`, `get_weather`

## Examples

``` r
if (FALSE) { # \dontrun{
sites <- data.frame(env = paste0("E", 1:200),
                    lat = runif(200, -30, -5),
                    lon = runif(200, -60, -40))

wth <- get_weather_resumable(
  env.id = sites$env, lat = sites$lat, lon = sites$lon,
  start.day = "2023-10-01", end.day = "2024-02-01",
  run.id = "trial2024")

## --- session dies at environment 137 ---
## In a NEW session, no arguments needed:
wth <- restart_from_log(run.id = "trial2024")

read_progress_log(run.id = "trial2024")
attr(wth, "failed")

## Straight into the characterization pipeline
wth <- processWTH(wth)
W   <- W_matrix(wth, env.id = "env", var.id = c("T2M", "PRECTOT"))
} # }
```
