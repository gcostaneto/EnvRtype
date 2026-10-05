# Kernel Models for Predicting Phenotypes across Multi-Environment Conditions

Fits Bayesian linear mixed models for multi-environment trials (MET)
using the kernels produced by
[`get_kernel()`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md).
The model is \$\$y = 1\mu + X\beta + \sum_j u_j + e\$\$ with \\u_j \sim
N(0, K_j \sigma^2_j)\\ and \\e \sim N(0, I\sigma^2_e)\\. Sampling runs
in the spectral basis of each kernel, so every kernel is used through
its eigendecomposition rather than as an n x n matrix.

The function returns genotype predictions, per-kernel variance
components, broad- and narrow-sense heritabilities, and – with
`keep_effects = TRUE` – posterior effect chains and breeder-facing
report tables (per-environment rankings, stability, selection
candidates).

Set `verbose = TRUE` to trace the six execution steps with per-step and
total timings; `verbose = FALSE` (default) is entirely silent.

## Usage

``` r
kernel_model(
  y,
  data = NULL,
  random = NULL,
  fixed = NULL,
  env,
  gid,
  verbose = FALSE,
  iterations = 1000,
  burnin = 200,
  thining = 10,
  tol = 1e-10,
  R2 = 0.5,
  digits = 4,
  seed = NULL,
  probs = c(0.025, 0.975),
  variance = c("realized", "scale"),
  scale_kernels = TRUE,
  keep_effects = FALSE,
  compute_CI = TRUE,
  R2_split = c("per_kernel", "shared"),
  prior_scale = c("eigen", "diag"),
  varcomp_method = c("postmean", "gibbs")
)
```

## Arguments

- y:

  character. Name of the phenotype column in `data`.

- data:

  data.frame with the environment, genotype and phenotype columns. Rows
  must be grouped by environment if any kernel is block-diagonal.

- random:

  list of kernels from
  [`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md).
  Each element is a list carrying `$decomp` (the spectral decomposition)
  and `$Type` (`"D"` or `"BD"`); `$Kernel` is optional and is not
  required by the sampler. Factored kernels (no `$Kernel`) and dense
  kernels may be mixed in one list. Names should follow the `KG_` /
  `KE_` / `KGE_` convention, which drives the variance-component roles.

- fixed:

  matrix (`n x p`) of fixed effects, or `NULL`.

- env, gid:

  character. Environment and genotype column names in `data`.

- verbose:

  logical or 0/1. `TRUE` traces the six execution steps with per-step
  and total timings; `FALSE` (default) is entirely silent – no
  [`cat()`](https://rdrr.io/r/base/cat.html), no
  [`message()`](https://rdrr.io/r/base/message.html). Results are
  identical either way.

- iterations:

  integer. Total Gibbs iterations. Default `1000`.

- burnin:

  integer. Iterations discarded before collecting draws. Default `200`.

- thining:

  integer. Keep every `thining`-th post-burn-in draw. Default `10`.

- tol:

  numeric. Eigenvalue cutoff. Must match `get_kernel(tol =)` or the
  cached decomposition is rejected and recomputed. Default `1e-10`.

- R2:

  numeric in (0, 1). Prior proportion of phenotypic variance explained
  by the random terms. Default `0.5`.

- digits:

  integer. Rounding for the variance-component tables.

- seed:

  integer or `NULL`. Sets the RNG state for reproducibility.

- probs:

  length-2 numeric. Quantiles for posterior prediction intervals.
  Default `c(0.025, 0.975)`.

- variance:

  `"realized"` (default) or `"scale"`. See *Details*.

- scale_kernels:

  logical. Normalise each kernel to `mean(diag(K)) = 1`. Must match
  `get_kernel(scale_kernels =)` for the cache to be reused. Default
  `TRUE`.

- keep_effects:

  logical. Retain per-draw effect chains. Required for
  `varcomp_method = "postmean"` and for the breeder-facing report tables
  (`$predictions`, `$genotype_effects`, `$env_summary`). Costs O(k \*
  draws \* n) memory. Default `FALSE`.

- compute_CI:

  logical. Compute posterior prediction intervals for `yHat`. Requires
  storing effect chains. Default `TRUE`.

- R2_split:

  `"per_kernel"` (default) or `"shared"`. Controls how the prior
  variance `R2` is distributed across kernels. See *Prior
  parameterisation*.

- prior_scale:

  `"eigen"` (default) or `"diag"`. Whether the prior scale is divided by
  the mean eigenvalue of each kernel, making kernels of differing rank
  comparable. See *Prior parameterisation*.

- varcomp_method:

  `"postmean"` (default) or `"gibbs"`. The estimator used for
  `$varcomp`. `"postmean"` needs `keep_effects = TRUE`; when unavailable
  it falls back to `"gibbs"` (reported in the trace when
  `verbose = TRUE`).

## Value

An object of class `"kernel_model"`, a list with:

- `yHat`:

  numeric vector of fitted/predicted values, length n.

- `yHat.CI`:

  matrix of posterior prediction intervals (`NULL` unless
  `compute_CI = TRUE`).

- `varE`:

  posterior mean residual variance.

- `random`:

  per-kernel posterior effects and variances.

- `fit`:

  the full sampler object: chains, kernel names, MCMC settings and –
  with `keep_effects = TRUE` – `Uchain`. (Named `$BGGE` in versions
  before v1.2.3.)

- `meta`:

  the aligned env / gid / y data frame.

- `VarComp`:

  variance-component table with confidence intervals.

- `varcomp`:

  validated components and heritabilities, class `"km_varcomp"`; see
  `varcomp_summary`.

- `varcomp_gibbs`:

  the same computed by the Gibbs-chain estimator, for comparison.

- `runtime`:

  sampler wall time in seconds.

With `keep_effects = TRUE` it additionally carries `predictions`,
`genetic_cor_env`, `genotype_effects`, `env_summary` and
`variance_proportions`.

## Details

Fits \$\$y = 1\mu + X\beta + \sum_j u_j + e\$\$ with \\u_j \sim N(0, K_j
\sigma^2\_{u_j})\\ and \\e \sim N(0, I\sigma^2_e)\\, by Gibbs sampling
in the spectral basis of each kernel. Writing \\K_j = U_j S_j U_j'\\,
the effects are sampled as \\u_j = U_j b_j\\ with \\b_j \sim N(0, S_j
\sigma^2\_{u_j})\\, so the cost per iteration scales with the retained
**rank** of each kernel, not with \\n^2\\. Trimming the eigenvalue tail
via `get_kernel(keep_var =)` is therefore the most direct way to speed
up sampling.

`variance = "scale"` reports the raw kernel-scale parameter `sigb[j]`
(original BGGE behaviour, inflated for low-rank incidence kernels such
as `KE_`). `variance = "realized"` reports `var(u_j)` averaged over
draws – what most users mean by a variance component. The raw scale is
always returned as `$sigb`.

## Original Version

Costa-Neto et al (2021)

EnvRtype v.1.2.3, Sep 2026

## Performance changes vs. the original

- **Cached eigendecomposition.** If `random` carries `$decomp` (attached
  by `get_kernel(decompose = TRUE)`) and the cache matches `tol` /
  `scale_kernels` / `n`, [`eigen()`](https://rdrr.io/r/base/eigen.html)
  is skipped entirely. Cross-validation folds and multiple chains reuse
  one decomposition.

- **No stored transposes.** The original kept both `U` and `t(U)` per
  kernel (doubling eigenvector memory) and computed `crossprod(tU, b)`.
  Now only `U` is stored and the back-transform is `U %*% b`.

- **Optional chain storage.** Full per-draw effect matrices
  (`k x nCum x n` doubles) are allocated only when needed
  (`keep_effects = TRUE` or `yHat.CI`). Otherwise running sums reduce
  memory from O(k\*nCum\*n) to O(k\*n).

- **Vectorised post-processing.** Row-wise
  [`apply()`](https://rdrr.io/r/base/apply.html) over draw matrices
  replaced by `rowSums`/`rowMeans` identities;
  [`sweep()`](https://rdrr.io/r/base/sweep.html) replaced by recycled
  arithmetic.

- **Scalar chains in a matrix** rather than a nested list, removing
  per-iteration copy-on-modify of the chain container.

## Prior parameterisation

Two defects in the original parameterisation are addressed by arguments:

**`R2_split`**. The original gave *every* kernel the full `R2` share, so
with `nk` kernels the priors collectively assert `nk * R2` of the
phenotypic variance – with `RNMDs` (4 kernels, `R2 = 0.5`) that is twice
`var(y)`. `"per_kernel"` (default) uses `R2 / nk` so total prior mass
stays at `R2` regardless of model size, keeping model comparisons on
equal footing. `"shared"` restores the original behaviour.

**`prior_scale`**. Kernels normalised to `mean(diag(K)) = 1` are equal
in scale but not in effective dimension: a rank-5 and a rank-750 kernel
with the same mean diagonal carry very different total variance.
`"eigen"` (default) divides the prior scale by `mean(s)` (= trace /
rank), making priors comparable across kernels of differing rank.

## References

Costa-Neto G, Galli G, Carvalho HF, Crossa J, Fritsche-Neto R (2021).
EnvRtype: a software to interplay enviromics and quantitative genomics
in agriculture. *G3*, 11(4), jkab040.

Granato I, Cuevas J, Luna-Vazquez F, Crossa J, Montesinos-Lopez O,
Burgueno J, Fritsche-Neto R (2018). BGGE: a new package for
genomic-enabled prediction incorporating genotype x environment
interaction models. *G3*, 8(9), 3039-3047.

## See also

[`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md),
`varcomp_summary`, `kernel_model_mc`, `kernel_cv`

## Author

Germano Costa Neto. Gibbs sampler adapted from Granato et al. (2018),
BGGE package. Refactored version.

## Examples

``` r
if (FALSE) { # \dontrun{
data("maizeYield"); data("maizeG"); data("maizeWTH")

ECs <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
                var.id = c("FRUE", "PETP", "SRAD", "T2M_MAX"),
                statistic = "mean")
K <- get_kernel(K_G = list(G = maizeG),
                K_E = list(W = env_kernel(env.data = ECs)[[2]]),
                data = maizeYield, model = "RNMM",
                env = "env", gid = "gid", y = "value")

## silent fit
fit <- kernel_model(y = "value", data = maizeYield, random = K,
                    env = "env", gid = "gid",
                    iterations = 5000, burnin = 1000)
fit
fit$varcomp

## traced fit, with effect chains for the report tables
fit2 <- kernel_model(y = "value", data = maizeYield, random = K,
                     env = "env", gid = "gid",
                     iterations = 5000, burnin = 1000,
                     keep_effects = TRUE, verbose = TRUE)
head(fit2$genotype_effects)

## leave-one-environment-out: kernels do not depend on y, so build once
for (e in unique(maizeYield$env)) {
  tr <- maizeYield; tr$value[tr$env == e] <- NA
  f <- kernel_model(y = "value", data = tr, random = K,
                    env = "env", gid = "gid", iterations = 2000, burnin = 500)
}
} # }
```
