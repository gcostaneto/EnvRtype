# Simulate a q x q Genetic Correlation Among Environments

Accessory to
[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md)
and [`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md).
Builds a valid \\q \times q\\ correlation matrix among environments, so
the same `C` can seed both simulators explicitly.

## Usage

``` r
sim_met_C(
  q,
  min.cor = 0,
  max.cor = 0.8,
  env_names = NULL,
  max.tries = 200L,
  seed = NULL
)
```

## Arguments

- q:

  integer. Number of environments (\>= 3).

- min.cor, max.cor:

  numeric. Range for off-diagonal correlations.

- env_names:

  character. Optional environment names.

- max.tries:

  integer. Rejection-sampling attempts before repair.

- seed:

  integer. RNG seed.

## Value

A \\q \times q\\ correlation matrix of class `"C_env"` with a logical
attribute `"repaired"` and an attribute `"requested"` recording the
target range.

## Details

Sampling off-diagonals uniformly from `[min.cor, max.cor]` does not
generally give a positive semi-definite matrix. The equicorrelated bound
is \\\rho \ge -1/(q-1)\\, so for \\q = 6\\ nothing below \\-0.2\\ is
attainable. A naive eigenvalue repair silently returns something else,
so this function errors up front when the request is infeasible and
rejection-samples otherwise, only falling back to eigenvalue clamping
(with a warning) if sampling fails.

## See also

[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md),
[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md),
[`env_cor`](https://gcostaneto.github.io/EnvRtype/reference/env_cor.md)

## Examples

``` r
C <- sim_met_C(q = 8, min.cor = 0.2, max.cor = 0.8, seed = 1)
range(C[lower.tri(C)])
#> [1] 0.2505481 0.7763708
```
