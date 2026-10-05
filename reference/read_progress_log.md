# Inspect the Progress Log of a Checkpointed Run

Returns the per-environment retrieval history written by
`get_weather_resumable` or `get_soil_resumable`: what succeeded, what
failed, when, and with which error message.

## Usage

``` r
read_progress_log(run.id, dir.path = NULL, all = FALSE, verbose = TRUE)
```

## Arguments

- run.id:

  character. Run directory name.

- dir.path:

  character. Parent directory. Default
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- all:

  logical. If `TRUE`, return every record including superseded retries.
  If `FALSE` (default), the latest record per environment.

- verbose:

  boolean. If `TRUE` (default) prints a progress banner.

## Value

A `data.frame` with `env`, `key`, `status` (`"ok"` or `"failed"`),
`time` and `message`.

## See also

`restart_from_log`

## Examples

``` r
if (FALSE) { # \dontrun{
lg <- read_progress_log(run.id = "zoning_grid")
table(lg$status)

## Why did they fail?
subset(lg, status == "failed")[, c("env", "message")]

## Full history including retries
read_progress_log(run.id = "zoning_grid", all = TRUE)
} # }
```
