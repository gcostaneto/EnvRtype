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

# S3 method for class 'sim_W'
env_cor(x, type = c("target", "realised"), ...)

# Default S3 method
env_cor(x, type = c("target", "realised"), ...)

# S3 method for class 'sim_envirome'
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
# \donttest{
met <- sim_met(maizeG, n_env = 6, seed = 1)
#> <sim_met>
#>   mechanism .............. kronecker
#>   lines x environments ... 150 x 6   (K is 150 x 150, C_env is 6 x 6)
#>   observations ........... 900  (n_rep = 1)
#>   h2 target .............. 0.312 - 0.696
#>   h2 realised (plot) ..... 0.297 - 0.557   max err 0.263
#>   h2 realised (line-mean). 0.297 - 0.557
#>   gcor target ............ 0.011 - 0.794
#>   gcor realised (whitened) 0.011 - 0.794   max err 0.00525
#>   exact .................. TRUE
env_cor(met, "target")
#>           E1        E2         E3        E4        E5         E6
#> E1 1.0000000 0.3981594 0.57409481 0.7935249 0.3040281 0.62195618
#> E2 0.3981594 1.0000000 0.74776418 0.1697140 0.5213390 0.10044408
#> E3 0.5740948 0.7477642 1.00000000 0.2137765 0.3088913 0.01071227
#> E4 0.7935249 0.1697140 0.21377653 1.0000000 0.3059104 0.69575268
#> E5 0.3040281 0.5213390 0.30889127 0.3059104 1.0000000 0.27227920
#> E6 0.6219562 0.1004441 0.01071227 0.6957527 0.2722792 1.00000000
#> attr(,"tries")
#> [1] 2
#> attr(,"repaired")
#> [1] FALSE
#> attr(,"requested")
#> min.cor max.cor 
#>     0.0     0.8 
#> attr(,"structure")
#> [1] "random"
#> attr(,"realised")
#>           E1        E2         E3        E4        E5         E6
#> E1 1.0000000 0.3990479 0.57413944 0.7935246 0.3044541 0.62199307
#> E2 0.3990479 1.0000000 0.74860862 0.1700671 0.5265888 0.10123990
#> E3 0.5741394 0.7486086 1.00000000 0.2137838 0.3100444 0.01083931
#> E4 0.7935246 0.1700671 0.21378385 1.0000000 0.3062818 0.69578108
#> E5 0.3044541 0.5265888 0.31004442 0.3062818 1.0000000 0.27224457
#> E6 0.6219931 0.1012399 0.01083931 0.6957811 0.2722446 1.00000000
#> attr(,"class")
#> [1] "C_env"  "matrix" "array" 
env_cor(met, "realised")
#>           E1        E2         E3        E4        E5         E6
#> E1 1.0000000 0.3990479 0.57413944 0.7935246 0.3044541 0.62199307
#> E2 0.3990479 1.0000000 0.74860862 0.1700671 0.5265888 0.10123990
#> E3 0.5741394 0.7486086 1.00000000 0.2137838 0.3100444 0.01083931
#> E4 0.7935246 0.1700671 0.21378385 1.0000000 0.3062818 0.69578108
#> E5 0.3044541 0.5265888 0.31004442 0.3062818 1.0000000 0.27224457
#> E6 0.6219931 0.1012399 0.01083931 0.6957811 0.2722446 1.00000000
# }
```
