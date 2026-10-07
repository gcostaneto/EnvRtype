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
  cal.nrep = 15L,
  marginal = c("gaussian", "right-skewed", "heavy-tailed", "bounded"),
  kernel = c("linear", "rbf"),
  bandwidth = NULL,
  na.frac = 0,
  family.sizes = NULL,
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

  numeric in \[0, 1). Share of columns collapsed onto drivers (also the
  within-family strength when `family.sizes` is used).

- n_blocks:

  integer. Collinear columns grouped into this many families.

- calibrate:

  logical. Solve for the internal variance share \\a\\ numerically so
  `noise` is on the outcome scale. Default `TRUE`.

- cal.nrep:

  integer. Monte-Carlo replicates used during calibration.

- marginal:

  character. Marginal distribution of the covariables: `"gaussian"`
  (default), `"right-skewed"`, `"heavy-tailed"` or `"bounded"`. Applied
  rank-preserving, so recovery is barely changed.

- kernel:

  character. Environmental kernel: `"linear"` (default, \\WW^\top\\) or
  `"rbf"` (Gaussian).

- bandwidth:

  numeric or `NULL`. RBF bandwidth; `NULL` uses the median heuristic.

- na.frac:

  numeric in \[0, 1). Fraction of entries set to `NA` in the returned
  \\W\\ (truth is computed on the complete matrix first).

- family.sizes:

  integer vector (optionally named). Column counts per covariable
  family; must sum to `n_var`. Overrides `n_blocks` redundancy and names
  the columns by family.

- var_names, env_names:

  character. Optional names for columns / rows.

- seed:

  integer. RNG seed. The caller's RNG stream is restored on exit.

- verbose:

  logical. Print a summary. Default `TRUE`.

## Value

A \\q \times n\_{var}\\ matrix of class `"sim_W"`. Attributes: `"C_env"`
(the target), `"K_W"` (the realised kernel) and `"explained"` (a list
with `r`, `r2`, per-PC variance shares, effective rank `eff_rank`, the
internal share `a_internal`, the `marginal`/`kernel`/`na.frac` used, and
the requested/target values).

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

**Realism knobs.** `marginal` maps each covariable to a non-Gaussian
shape (rank-preserving), `kernel = "rbf"` builds a Gaussian instead of a
linear kernel, `family.sizes` groups columns into named covariable
families (as in real enviromes), and `na.frac` injects missing values
into the returned \\W\\ *after* the truth is computed, so downstream
cleaning (e.g.
[`W_matrix`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md)
imputation) can be tested against a known answer.

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

## 1. An envirome that recovers ~70% of C (noise is on the outcome scale)
W <- sim_W(C, noise = 0.3, n_var = 50, seed = 1)
attr(W, "explained")$r

## 2. Perfect vs pure-noise enviromes
W_good  <- sim_W(C, noise = 0.0, n_var = 50)
W_noise <- sim_W(C, noise = 1.0, n_var = 50)

## 3. Redundant envirome: many collinear columns, low effective rank
W <- sim_W(C, noise = 0.3, n_var = 200, collinearity = 0.8, n_blocks = 5)

## 4. Non-Gaussian covariables (rank-preserving marginal transform)
W <- sim_W(C, noise = 0.3, n_var = 60, marginal = "right-skewed")
W <- sim_W(C, noise = 0.3, n_var = 60, marginal = "heavy-tailed")
W <- sim_W(C, noise = 0.3, n_var = 60, marginal = "bounded")

## 5. Gaussian RBF kernel instead of the linear WW'
W <- sim_W(C, noise = 0.3, n_var = 60, kernel = "rbf")
W <- sim_W(C, noise = 0.3, n_var = 60, kernel = "rbf", bandwidth = 2)

## 6. Inject missing values AFTER the truth is computed
W <- sim_W(C, noise = 0.3, n_var = 60, na.frac = 0.1)
sum(is.na(W))

## 7. Named covariable families (sizes must sum to n_var)
W <- sim_W(C, n_var = 12, family.sizes = c(temp = 4, rain = 5, rad = 3))
colnames(W)

## 8. Custom names and a faster/looser calibration
W <- sim_W(C, noise = 0.4, n_var = 30, cal.nrep = 5,
           env_names = paste0("Loc", 1:10),
           var_names = paste0("bio", 1:30), seed = 1)

## 9. Raw (uncalibrated) internal variance share
W <- sim_W(C, noise = 0.5, n_var = 50, calibrate = FALSE)
plot(W)
} # }
```
