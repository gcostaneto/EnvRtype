# Simulate a Multi-Environment Trial with a Known Genetic Architecture

Core simulator. Generates genotype-by-environment phenotypes whose truth
is known exactly: a genomic kinship \\K\\ among lines, a genetic
(type-B) correlation \\C\\ among environments, and a per-environment
heritability. It is the reference-answer generator against which
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
and
[`scan_untested_envs`](https://gcostaneto.github.io/EnvRtype/reference/scan_untested_envs.md)
can be checked, and it pairs with
[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md) to
study how well an envirome recovers \\C\\.

## Usage

``` r
sim_met(
  K = NULL,
  n_env = NULL,
  min.h2 = 0.2,
  max.h2 = 0.8,
  min.cor = 0,
  max.cor = 0.8,
  C = NULL,
  h2 = NULL,
  n_rep = 1,
  sg2 = 1,
  exact = TRUE,
  env_effects = NULL,
  subset = NULL,
  subset.by = c("cell", "gid_env"),
  K.weights = NULL,
  W = NULL,
  reaction.norm = !is.null(W),
  rn.main = 0.5,
  gxe.specific = 0,
  markers = NULL,
  n.qtl = NULL,
  rep.var = 0,
  spatial = FALSE,
  field.dim = NULL,
  spatial.rho = 0.6,
  spatial.var = 1,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- K:

  n x n genomic kinship AMONG LINES (e.g. `maizeG`), or a named *list*
  of such matrices for multiple variance components (combined with
  `K.weights`). Made positive definite internally, so a singular
  marker-based kinship is fine. May be `NULL` when `markers` is given.

- n_env:

  integer. Number of environments; inferred from `C` if given.

- min.h2, max.h2:

  numeric. Plot-basis heritability range, used when `h2` is not
  supplied. Must satisfy \\0 \< min.h2 \le max.h2 \< 1\\.

- min.cor, max.cor:

  numeric. Genetic correlation range among environments, used when `C`
  is not supplied. Mutually exclusive with `C`.

- C:

  q x q correlation AMONG ENVIRONMENTS (see
  [`sim_met_C`](https://gcostaneto.github.io/EnvRtype/reference/sim_met_C.md)).
  Mutually exclusive with `min.cor`/`max.cor`.

- h2:

  numeric vector of per-environment heritabilities (length 1 or
  `n_env`). Mutually exclusive with `min.h2`/`max.h2`.

- n_rep:

  integer. Replicates per genotype x environment cell.

- sg2:

  numeric. Genetic variance per environment (length 1 or `n_env`).

- exact:

  logical. Orthonormalise the draw so realised moments match the
  targets. Default `TRUE`.

- env_effects:

  numeric. Optional fixed environment means \\\mu_j\\.

- subset:

  numeric in (0, 1\]. Fraction of cells retained, to induce
  unbalancedness.

- subset.by:

  "cell" or "gid_env". Whether missingness is drawn per observation or
  per genotype-by-environment cell.

- K.weights:

  numeric. Non-negative weights (rescaled to sum 1) for a list of
  kinships in `K`; ignored for a single matrix. Default equal weights.

- W:

  q x k envirome matrix. When supplied, genetic values are generated as
  reaction norms and the environment correlation is induced from `W`.

- reaction.norm:

  logical. Use the reaction-norm mechanism. Defaults to `TRUE` when `W`
  is supplied.

- rn.main:

  numeric in \[0, 1\]. Share of genetic variance assigned to the main
  (across-environment) effect under the reaction-norm mechanism; the
  rest is reaction-norm interaction. Default 0.5.

- gxe.specific:

  numeric in \[0, 1). Share of genetic variance made line-independent
  and environment-specific, breaking separability. Default 0.

- markers:

  n x m marker matrix. When supplied, a marker/QTL mechanism replaces
  the infinitesimal draw and `K` is derived from the markers if not
  given.

- n.qtl:

  integer. Number of markers acting as QTL (sampled) when `markers` is
  used; `NULL` uses all markers.

- rep.var:

  numeric. Variance of a random replicate/block effect added to
  phenotypes. Default 0.

- spatial:

  logical. Add an AR1\\\times\\AR1 spatial field to each environment.
  Default `FALSE`.

- field.dim:

  integer length-2 `c(nrow, ncol)`. Field layout per environment;
  defaults to a near-square grid sized to the plots.

- spatial.rho:

  numeric. Lag-1 spatial autocorrelation. Default 0.6.

- spatial.var:

  numeric. Variance of the spatial field. Default 1.

- seed:

  integer. RNG seed. The caller's RNG stream is restored on exit.

- verbose:

  logical. Print a summary. Default `TRUE`.

## Value

An object of class `"sim_met"`: a list with

- `data`:

  long `data.frame` with `env`, `gid`, `rep`, `value`.

- `C_env`:

  the target environment correlation, with the realised correlation
  attached as attribute `"realised"`.

- `truth`:

  list of everything known: genetic values `g`, target/realised
  correlations, target/realised heritabilities on both bases, variances,
  means.

- `K`:

  the positive-definite kinship actually used.

## Details

**Model.** For \\n\\ lines and \\q\\ environments, genetic values are
drawn from a separable (Kronecker) Gaussian \$\$\mathrm{vec}(G) \sim
\mathcal{N}\\\left(0,\\ K \otimes \Sigma_g\right), \qquad \Sigma_g =
D^{1/2}\\ C\\ D^{1/2},\\ D = \mathrm{diag}(\sigma^2_g),\$\$ realised as
\\G = L_K\\ U\\ L_C^{\top}\\ with \\L_K L_K^{\top}=K\\, \\L_C
L_C^{\top}=\Sigma_g\\. When `exact = TRUE` the driver \\U\\ is
orthonormalised so that \\U^{\top}U = nI\\, making the realised moments
match their targets instead of merely scattering around them.

**Phenotypes.** For replicate \\r\\, \$\$y\_{ijr} = \mu_j + g\_{ij} +
\varepsilon\_{ijr}, \qquad
\varepsilon\_{ijr}\sim\mathcal{N}(0,\sigma^2\_{e,j}),\quad
\sigma^2\_{e,j} = \sigma^2\_{g,j}\\\frac{1-h^2_j}{h^2_j}.\$\$

**Heritability, two bases.** The plot-basis value \\h^2 =
V_g/(V_g+V_e)\\ is what `min.h2`/`max.h2` target; the line-mean basis
\\h^2/(h^2 + (1-h^2)/n\_{rep})\\ is also reported.

**Realised environment correlation is measured on the K-whitened
scale.** Because lines are related through \\K\\, `cor(G)` across lines
is not \\C\\; forcing it to be is wrong. The realised correlation is
defined as \\\mathrm{cor}(L_K^{-1} G)\\ and stored for validation.

**The two matrices are different objects.** `K` is \\n\times n\\ among
lines; the environment correlation is \\q\times q\\. A genomic kinship
passed where an environment correlation is expected is the most common
error and is rejected by an elementwise unit-diagonal check.

**Alternative genetic mechanisms.** Besides the default separable draw,
three opt-in mechanisms generate the genetic values:

- *Reaction norm* (supply `W`): genetic values are built causally from
  an envirome, \\g\_{ij} = m_i + \sum_k w\_{jk} b\_{ik}\\ with a main
  effect \\m\\ and line sensitivities \\b\_{\cdot k}\\ both
  \\K\\-structured. The environment correlation is then *induced* by `W`
  (and reported), not imposed. `rn.main` splits variance between the
  main effect and the reaction-norm interaction.

- *Marker-based* (supply `markers`): QTL effects are sampled per
  environment with the target correlation `C`, giving large-effect
  architectures instead of the infinitesimal draw.

- *Multiple kinships* (supply a *list* for `K` with `K.weights`):
  independent additive, dominance or epistatic components are summed.

`gxe.specific` adds a line-independent, environment-specific deviation
that makes the covariance *non-separable*. `rep.var` and `spatial` (an
AR1\\\times\\AR1 field) add trial realism to the phenotypes.

## References

Costa-Neto, G., et al. (2023). Envirome-wide associations enhance
multi-environment prediction. *G3* 13(2), jkac313.

## See also

[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md),
[`sim_met_C`](https://gcostaneto.github.io/EnvRtype/reference/sim_met_C.md),
[`env_cor`](https://gcostaneto.github.io/EnvRtype/reference/env_cor.md),
[`sim_markers`](https://gcostaneto.github.io/EnvRtype/reference/sim_markers.md),
[`sim_envirome`](https://gcostaneto.github.io/EnvRtype/reference/sim_envirome.md),
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md),
[`scan_untested_envs`](https://gcostaneto.github.io/EnvRtype/reference/scan_untested_envs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(maizeG)

## 1. Baseline: random C and random h2 within a range
met <- sim_met(maizeG, n_env = 10, min.cor = 0.2, max.cor = 0.8,
               min.h2 = 0.3, max.h2 = 0.7, seed = 1)
head(met$data)
env_cor(met, "realised")

## 2. Supply an explicit environment correlation (e.g. a year-like series)
C   <- sim_met_C(q = 8, structure = "ar1", rho = 0.6)
met <- sim_met(maizeG, C = C, h2 = 0.5, seed = 1)

## 3. Fixed per-environment heritabilities and environment means
met <- sim_met(maizeG, n_env = 4, h2 = c(0.3, 0.5, 0.6, 0.8),
               env_effects = c(2, 4, 6, 8), seed = 1)

## 4. Replicates and an unbalanced design (retain 70% of g x e cells)
met <- sim_met(maizeG, n_env = 6, n_rep = 3,
               subset = 0.7, subset.by = "gid_env", seed = 1)

## 5. Several genetic variance components (additive + dominance + epistasis)
met <- sim_met(K = list(add = maizeG, dom = maizeG, epi = maizeG),
               K.weights = c(0.6, 0.3, 0.1), n_env = 8, seed = 1)

## 6. Causal reaction-norm: C is INDUCED by the envirome, not imposed
Wenv <- scale(matrix(rnorm(10 * 6), 10, 6))
rn   <- sim_met(maizeG, W = Wenv, rn.main = 0.4, seed = 1)
env_cor(rn, "target")

## 7. Marker / QTL architecture instead of the infinitesimal draw
X   <- matrix(rbinom(200 * 500, 2, 0.3), 200, 500)
qtl <- sim_met(markers = X, n_env = 6, n.qtl = 30, seed = 1)

## 8. Non-separable GxE: 25% of genetic variance is environment-specific
met <- sim_met(maizeG, n_env = 6, gxe.specific = 0.25, seed = 1)

## 9. Field realism: replicate/block noise plus an AR1 x AR1 spatial trend
met <- sim_met(maizeG, n_env = 4, n_rep = 2, rep.var = 0.5,
               spatial = TRUE, field.dim = c(20, 15),
               spatial.rho = 0.7, spatial.var = 1.5, seed = 1)

## 10. Pair with a simulated envirome and diagnose the recovery
W <- sim_W(met$C_env, noise = 0.3, n_var = 50, seed = 1)
plot(met)

## 11. End-to-end: simulate kinship + MET + envirome, then fit a model
##     Simulate a genomic kinship KG from markers (VanRaden)
X <- matrix(rbinom(150 * 1000, 2, 0.3), 150, 1000)
rownames(X) <- paste0("G", seq_len(150))        # gid names carry through
Z <- scale(X); Z[is.na(Z)] <- 0
Kg <- tcrossprod(Z) / ncol(Z) + diag(1e-6, 150) # positive-definite kinship
dimnames(Kg) <- list(rownames(X), rownames(X))
KG <- list(G = Kg)

##     Simulate the MET truth on that kinship (one replicate)
met <- sim_met(Kg, n_env = 10, n_rep = 1,
               min.cor = 0.3, max.cor = 0.8,
               min.h2 = 0.3, max.h2 = 0.7, seed = 1)
df  <- as.data.frame(met$data)                  # base data.frame for get_kernel()

##     Simulate a redundant, noisier envirome on the same environments
W2  <- sim_W(met$C_env, noise = 0.6, n_var = 200,
             collinearity = 0.8, n_blocks = 5, seed = 1)
KE2 <- list(W = env_kernel(env.data = as.matrix(W2))[[2]])

##     Build reaction-norm kernels and fit the Bayesian model
K2  <- get_kernel(K_G = KG, K_E = KE2, data = df, model = "RNMM",
                  env = "env", gid = "gid", y = "value")
fit <- kernel_model(y = "value", data = df, random = K2,
                    env = "env", gid = "gid",
                    iterations = 5000, burnin = 1000)

## 12. Phenotypes for contrasting GENETIC diversity (envirome fixed)
##     Low diversity: few rare markers, strong structure; high: many common.
E    <- sim_envirome(12, 40, structure = "clusters", n_cluster = 3, seed = 1)
X_lo <- sim_markers(200, 500,  maf = "rare",   n_subpop = 5, Fst = 0.2, seed = 1)
X_hi <- sim_markers(200, 4000, maf = "common", n_subpop = 1, Fst = 0.0, seed = 1)
met_lo <- sim_met(markers = X_lo, C = attr(E, "C_env"), n.qtl = 80, seed = 1)
met_hi <- sim_met(markers = X_hi, C = attr(E, "C_env"), n.qtl = 80, seed = 1)
c(low = var(met_lo$data$value), high = var(met_hi$data$value))

## 13. Phenotypes for contrasting ENVIROME diversity (genetics fixed)
##     Diverse environments induce lower cross-environment correlation (more GxE).
X    <- sim_markers(200, 3000, n_subpop = 3, Fst = 0.08, seed = 2)
E_lo <- sim_envirome(12, 40, structure = "gradient", redundancy = 0.7,
                     n_blocks = 2, seed = 2)
E_hi <- sim_envirome(12, 40, structure = "random", redundancy = 0, seed = 2)
met_eLo <- sim_met(markers = X, C = attr(E_lo, "C_env"), n.qtl = 80, seed = 2)
met_eHi <- sim_met(markers = X, C = attr(E_hi, "C_env"), n.qtl = 80, seed = 2)
cor(met_eLo$truth$g)[1, 2]   # high  cross-environment genetic correlation
cor(met_eHi$truth$g)[1, 2]   # low   correlation -> stronger GxE

## 14. Full 2 x 2 factorial of genetic x envirome diversity
scenario <- function(gen, env, seed = 1) {
  X <- if (gen == "low")
         sim_markers(200, 500, maf = "rare", n_subpop = 5, Fst = 0.2,
                     seed = seed, verbose = FALSE)
       else
         sim_markers(200, 4000, maf = "common", n_subpop = 1, Fst = 0,
                     seed = seed, verbose = FALSE)
  E <- if (env == "low")
         sim_envirome(12, 40, structure = "gradient", redundancy = 0.7,
                      n_blocks = 2, seed = seed, verbose = FALSE)
       else
         sim_envirome(12, 40, structure = "random", redundancy = 0,
                      seed = seed, verbose = FALSE)
  met <- sim_met(markers = X, C = attr(E, "C_env"),
                 n.qtl = 80, min.h2 = 0.3, max.h2 = 0.6, seed = seed)
  Cg  <- cor(met$truth$g)
  data.frame(gen = gen, env = env,
             gxe_cor   = mean(Cg[lower.tri(Cg)]),
             pheno_var = var(met$data$value))
}
do.call(rbind, list(scenario("low",  "low"),  scenario("low",  "high"),
                    scenario("high", "low"),  scenario("high", "high")))

## 15. Sweep population structure (Fst) and record phenotype variance
E <- sim_envirome(10, 40, structure = "clusters", n_cluster = 3, seed = 3)
do.call(rbind, lapply(c(0, 0.05, 0.1, 0.2), function(f) {
  X   <- sim_markers(200, 3000, n_subpop = 4, Fst = f, seed = 3, verbose = FALSE)
  met <- sim_met(markers = X, C = attr(E, "C_env"), n.qtl = 80, seed = 3)
  data.frame(Fst = f, eff_dim = attr(X, "diversity")["eff_dim"],
             pheno_var = var(met$data$value))
}))

## 16. Reaction-norm phenotypes driven directly by a diverse envirome
X    <- sim_markers(150, 2500, n_subpop = 2, Fst = 0.06, seed = 6)
E2   <- sim_envirome(10, 25, structure = "gradient",
                     marginal = "right-skewed", seed = 6)
met2 <- sim_met(markers = X, W = as.matrix(E2),
                reaction.norm = TRUE, rn.main = 0.6, seed = 6)
met2$truth$mechanism
} # }
```
