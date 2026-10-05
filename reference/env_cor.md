# Extract the Correlation Among Environments

Accessory generic that pulls the environment correlation out of a
[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md)
object or a bare matrix, so a hand-off to
[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md)
never depends on internal structure.

## Usage

``` r
env_cor(x, type = c("target", "realised"), ...)

# S3 method for class 'sim_met'
env_cor(x, type = c("target", "realised"), ...)

# S3 method for class 'matrix'
env_cor(x, type = c("target", "realised"), ...)

# Default S3 method
env_cor(x, type = c("target", "realised"), ...)
```

## Arguments

- x:

  an object (`sim_met`, matrix, ...).

- type:

  "target" (the specified correlation) or "realised" (what the
  simulation actually produced, on the K-whitened scale).

- ...:

  unused.

## Value

A \\q \times q\\ correlation matrix.

## See also

[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md),
[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md)

## Examples

``` r
if (FALSE) { # \dontrun{
met <- sim_met(maizeG, n_env = 6, seed = 1)
env_cor(met, "target")
env_cor(met, "realised")
} # }
```
