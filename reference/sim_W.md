# Simulate an Envirome with a Known Relationship to an Environment Correlation

Core simulator. Given a target correlation \\C\_{env}\\ among
environments, it builds a \\q \times n\_{var}\\ covariable matrix \\W\\
whose environmental kernel recovers a *known* share of \\C\_{env}\\. It
answers: if the true similarity among my environments is \\C\_{env}\\,
what would an envirome of this size and quality look like, and how much
of that similarity would it recover? Pairs with
[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md).

## Usage

``` r
sim_W(
  C_env,
  noise = 0,
  n_var = 50,
  collinearity = 0,
  n_blocks = 1,
  calibrate = TRUE,
  var_names = NULL,
  env_names = NULL,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- C_env:

  q x q correlation AMONG ENVIRONMENTS, from any source (e.g.
  [`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md)`$C_env`,
  [`env_cor`](https://gcostaneto.github.io/EnvRtype/reference/env_cor.md),
  `cov2cor(env_kernel(W)$envCov)`, or hand-specified). NOT a kinship.

- noise:

  numeric in \[0, 1\]. With `calibrate = TRUE`, the expected
  off-diagonal correlation between the simulated kernel and `C_env` is
  `1 - noise`.

- n_var:

  integer. Number of covariables (columns of \\W\\).

- collinearity:

  numeric in \[0, 1). Share of columns collapsed onto drivers.

- n_blocks:

  integer. Collinear columns grouped into this many families.

- calibrate:

  logical. Solve for the internal variance share \\a\\ numerically so
  `noise` is on the outcome scale. Default `TRUE`.

- var_names, env_names:

  character. Optional names for columns / rows.

- seed:

  integer. RNG seed.

- verbose:

  logical. Print a summary. Default `TRUE`.

## Value

A \\q \times n\_{var}\\ matrix of class `"sim_W"`. Attributes: `"C_env"`
(the target), `"K_W"` (the realised kernel) and `"explained"` (a list
with `r`, `r2`, per-PC variance shares, effective rank `eff_rank`, the
internal share `a_internal`, and the requested/target values).

## Details

**Construction.** With \\C\_{env}^{1/2}\\ the symmetric PSD square root
of the target, an orthonormal driver \\O\\ (\\q\times n\_{var}\\) and
iid noise \\E\\, \$\$W = \sqrt{1-a}\\ C\_{env}^{1/2} O \\ + \\
\sqrt{a}\\ E,\$\$ and the environmental kernel follows the
[`env_kernel`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md)
convention \$\$K_W = \frac{W W^{\top}}{\mathrm{tr}(W W^{\top})/q}.\$\$
Recovery is measured by the off-diagonal Pearson correlation \\r =
\mathrm{cor}\big(\mathrm{offdiag}(K_W),\\
\mathrm{offdiag}(C\_{env})\big)\\.

**Calibration.** When `calibrate = TRUE`, `noise` is on the *outcome*
scale: the internal variance share \\a\\ is solved by
[`uniroot`](https://rdrr.io/r/stats/uniroot.html) so that \\E\[r\] = 1 -
\mathrm{noise}\\. Because \\r\\ decreases in \\a\\, `noise = 0` gives
\\a = 0\\ (full recovery) and `noise = 1` gives \\a = 1\\ (pure noise).

**Collinearity.** With `collinearity > 0`, a share of columns is
collapsed onto `n_blocks` driver families, mimicking the redundancy of
real enviromes and lowering the effective rank.

**`C_env` is an environment correlation, not a kinship.** It must be
\\q\times q\\ with unit diagonal. Passing a genomic kinship is rejected.

## References

Costa-Neto, G., et al. (2023). Envirome-wide associations enhance
multi-environment prediction. *G3* 13(2), jkac313.

## See also

[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md),
[`sim_met_C`](https://gcostaneto.github.io/EnvRtype/reference/sim_met_C.md),
[`sim_W_grid`](https://gcostaneto.github.io/EnvRtype/reference/sim_W_grid.md),
[`env_cor`](https://gcostaneto.github.io/EnvRtype/reference/env_cor.md),
[`env_kernel`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
C <- sim_met_C(q = 10, min.cor = 0.2, max.cor = 0.8, seed = 1)

## An envirome that recovers ~70% of C
W <- sim_W(C, noise = 0.3, n_var = 50, seed = 1)
attr(W, "explained")$r

## Redundant envirome: many collinear columns, low effective rank
W2 <- sim_W(C, noise = 0.3, n_var = 200, collinearity = 0.8, n_blocks = 5)
} # }
```
