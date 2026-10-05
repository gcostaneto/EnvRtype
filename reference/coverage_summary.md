# Summarise coverage of a scanned region

Tabulates how much of a scanned area falls in each envelope class, and
which covariates most often drive sites outside. The second table
answers "what should I go and measure".

## Usage

``` r
coverage_summary(scan)
```

## Arguments

- scan:

  a \`scan_untested_envs\` / \[grid_scan()\] result.

## Value

list of \`by_position\` and \`by_covariate\` data.frames.
