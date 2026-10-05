# Resume an Interrupted get_weather() or get_soil() Run

Reads the manifest and progress log written by `get_weather_resumable`
or `get_soil_resumable`, skips every environment already retrieved,
fetches those still missing or failed, and returns the assembled result.

Safe to call repeatedly. If everything already succeeded it fetches
nothing and simply assembles what is on disk.

## Usage

``` r
restart_from_log(
  run.id,
  dir.path = NULL,
  retry.failed = TRUE,
  sleep = NULL,
  tries = NULL,
  verbose = TRUE
)
```

## Arguments

- run.id:

  character. The run directory name used originally.

- dir.path:

  character. Parent directory. Default
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- retry.failed:

  logical. Retry environments previously logged as failed. Default
  `TRUE`. Set `FALSE` to assemble only the successes – useful when
  failures are permanent, such as coordinates in the ocean.

- sleep:

  numeric or `NULL`. Override the stored delay. `NULL` (default) reuses
  the manifest value.

- tries:

  integer or `NULL`. Override the stored retry count.

- verbose:

  logical. Print progress. Default `TRUE`.

## Value

A `data.frame` as returned by the original function, with attributes
`run.dir`, `n_ok`, `n_failed` and `failed`.

## Details

All retrieval parameters come from `manifest.json`, so a resumed run
necessarily uses the same coordinates, dates and variables as the
original. You cannot accidentally resume with different settings – that
is what the stored fingerprint prevents.

**Permanently failing environments.** Some failures never resolve: a
point in the ocean has no SoilGrids profile. Repeated calls will retry
them every time. Once you have confirmed a failure is permanent, pass
`retry.failed = FALSE`.

## See also

`get_weather_resumable`, `get_soil_resumable`, `read_progress_log`

## Examples

``` r
if (FALSE) { # \dontrun{
## A long soil job dies at environment 137 of 200.
## In a brand-new session -- no need to remember the arguments:
soil <- restart_from_log(run.id = "zoning_grid")

## Inspect what failed and why
lg <- read_progress_log(run.id = "zoning_grid")
subset(lg, status == "failed")[, c("env", "message")]

## Give up on permanent failures
soil <- restart_from_log(run.id = "zoning_grid", retry.failed = FALSE)

## Weather resumes identically
wth <- restart_from_log(run.id = "trial2024")
} # }
```
