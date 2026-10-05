# Envirotype Frequency Matrix, Ready for the Kernel Pipeline

The envirotype counterpart of
[`W_matrix`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md).
Runs
[`env_typing`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md)
and returns, unconditionally, an environment x envirotype numeric matrix
suitable for
[`env_kernel`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md)
and
[`env_cluster`](https://gcostaneto.github.io/EnvRtype/reference/env_cluster.md)
– regardless of which of the four
[`env_typing()`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md)
branches ran. Where `W_matrix` summarises continuous covariables,
`T_matrix` summarises the *frequencies* of discrete environmental types
(envirotypes), the "enviromic assembly" building block of Costa-Neto,
Crossa & Fritsche-Neto (2021).

## Usage

``` r
T_matrix(
  env.data,
  var.id = NULL,
  env.id = NULL,
  cardinals = NULL,
  days.id = NULL,
  time.window = NULL,
  names.window = NULL,
  quantiles = NULL,
  id.names = NULL,
  by.interval = FALSE,
  joint = FALSE,
  joint.by.interval = TRUE,
  envirotype_mining = FALSE,
  k.range = 2:10,
  nstart = 25,
  iter.max = 100,
  seed = 1234,
  combine.features = FALSE,
  combine.max.order = 2L,
  ratio = TRUE,
  center = FALSE,
  scale = FALSE,
  sd.tol = 10,
  tol = 0.001,
  QC = FALSE,
  impute = c("none", "mean", "drop"),
  drop.constant = TRUE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame of processed weather data, or an object already returned by
  [`env_typing`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md)
  (any of its four shapes), in which case it is reshaped rather than
  recomputed.

- var.id, env.id, cardinals, days.id, time.window, names.window,
  quantiles:

  Passed to
  [`env_typing`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md).

- id.names, by.interval, joint, joint.by.interval:

  Passed to
  [`env_typing`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md).

- envirotype_mining, k.range, nstart, iter.max, seed:

  Passed to
  [`env_typing`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md).

- combine.features, combine.max.order:

  Passed to
  [`env_typing`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md).

- ratio:

  logical. Convert counts to within-environment relative frequencies.
  Default `TRUE` – see Details.

- center, scale:

  logical. Passed to the shared scaler. **Both default to `FALSE`**,
  unlike
  [`W_matrix`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md)
  – see Details.

- sd.tol, tol, QC, impute:

  Passed to the shared scaler.

- drop.constant:

  logical. Drop envirotype columns that are constant across environments
  (they contribute nothing to a kernel). Default `TRUE`.

- verbose:

  logical.

## Value

A numeric matrix of class `c("T_matrix", "matrix", "array")`, with
`rownames` = environments and `colnames` = envirotypes. Attributes:
`"typologies"` (the long frame), `"envirotype_description"`,
`"mining_summary"`, `"plot"`, `"ratio"`, `"dropped"`, and the
`.env_w_scale()` attributes when scaling is requested.

## Details

**Why `ratio = TRUE` by default.** Raw envirotype counts are days, so an
environment with a longer season has larger counts in every column. A
kernel built on raw counts is then partly a kernel on season length.
Row-normalising makes each environment a composition – the share of its
season spent in each envirotype – which is what "envirotype frequency"
is meant to convey.

**Why `center = scale = FALSE` by default.** This deliberately differs
from
[`W_matrix`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md),
where `center = scale = TRUE`. After `ratio = TRUE` the columns are
already on a common \[0, 1\] scale and each row sums to 1. Scaling a
composition to unit variance per column inflates rare envirotypes – a
type present in 2% of days in one environment and 0% elsewhere becomes a
high-leverage column. If you want the `W_matrix` convention, set them
explicitly; the choice is yours but it should be a choice.

**Compositional caveat.** With `ratio = TRUE` the rows sum to 1, so the
columns are linearly dependent and the matrix is rank-deficient by
exactly one. This is harmless for
[`env_kernel`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md)
(a Gaussian/linear kernel is still positive semi-definite) but will
produce one zero eigenvalue. Do not "fix" it by dropping a column – that
makes the result depend on which column you dropped.

## References

Costa-Neto, G., Crossa, J., & Fritsche-Neto, R. (2021). Enviromic
assembly increases accuracy and reduces costs of the genomic prediction
for yield plasticity in maize. *Frontiers in Plant Science* 12, 717552.

## See also

[`W_matrix`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md)
for the covariable counterpart,
[`env_typing`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md)
for the underlying typology,
[`env_kernel`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md)
and
[`env_cluster`](https://gcostaneto.github.io/EnvRtype/reference/env_cluster.md)
for the consumers.

## Examples

``` r
# \donttest{
data(maizeWTH)
Tm <- T_matrix(maizeWTH, var.id = c("T2M", "PRECTOT"), env.id = "env")
#> ---------------------------------------------------------------
#> T_matrix -- builds the envirotype (T) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - mining envirotypes with env_typing()
#>   - normalising envirotype counts to within-environment shares
#> <T_matrix>
#>   environments x envirotypes . 5 x 7
#>   row-normalised (ratio) ..... TRUE
#>   row sums ................... 1.0000 - 1.0000
dim(Tm)
#> [1] 5 7
K <- env_kernel(env.data = Tm, is.scaled = TRUE)
#> ---------------------------------------------------------------
#> env_kernel -- builds environmental relatedness kernels
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing a single environmental kernel
# }
```
