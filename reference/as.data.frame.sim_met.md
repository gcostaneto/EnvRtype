# Coerce a Simulated MET to a data.frame

Returns the long phenotype table, so a `sim_met` object can be handed
straight to modelling functions or inspected as a plain frame.

## Usage

``` r
# S3 method for class 'sim_met'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  a `sim_met` object.

- row.names, optional:

  passed to the default method for consistency; unused.

- ...:

  ignored.

## Value

A `data.frame` with one row per plot and columns `env`, `gid`, `rep` and
`value`.

## See also

[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md)
