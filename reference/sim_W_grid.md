# Sweep Envirome Size and Quality (Power Analysis)

Accessory to
[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md).
Repeatedly simulates enviromes across a grid of `noise`, `n_var` and
`collinearity`, returning the mean recovery and effective rank at each
design point.

## Usage

``` r
sim_W_grid(
  C_env,
  noise = seq(0, 1, 0.25),
  n_var = c(10, 50, 200),
  collinearity = c(0, 0.5, 0.9),
  n_rep = 10,
  calibrate = TRUE,
  cal.nrep = 15L,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- C_env:

  q x q correlation AMONG ENVIRONMENTS (see
  [`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md)).

- noise, n_var, collinearity:

  numeric vectors defining the grid.

- n_rep:

  integer. Replicate simulations per design point.

- calibrate:

  logical. Passed to
  [`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md).

- cal.nrep:

  integer. Calibration replicates per
  [`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md)
  call.

- seed:

  integer. RNG seed; each design point gets a reproducible sub-seed.

- verbose:

  logical. Print per-design-point progress. Default `TRUE`.

## Value

A `data.frame` with one row per design point: `noise`, `n_var`,
`collinearity`, mean recovery `r` and its SD `r_sd`, `r2`, and mean
`eff_rank`.

## See also

[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md),
[`sim_met_C`](https://gcostaneto.github.io/EnvRtype/reference/sim_met_C.md)

## Examples

``` r
if (FALSE) { # \dontrun{
C <- sim_met_C(q = 10, min.cor = 0.2, max.cor = 0.8, seed = 1)
grid <- sim_W_grid(C, noise = c(0, 0.3, 0.6), n_var = c(20, 100),
                   collinearity = c(0, 0.8), n_rep = 5)
} # }
```
