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
  K,
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
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- K:

  n x n genomic kinship AMONG LINES (e.g. `maizeG`). Made positive
  definite internally, so a singular marker-based kinship is fine.

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

- seed:

  integer. RNG seed.

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

## References

Costa-Neto, G., et al. (2023). Envirome-wide associations enhance
multi-environment prediction. *G3* 13(2), jkac313.

## See also

[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md),
[`sim_met_C`](https://gcostaneto.github.io/EnvRtype/reference/sim_met_C.md),
[`env_cor`](https://gcostaneto.github.io/EnvRtype/reference/env_cor.md),
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md),
[`scan_untested_envs`](https://gcostaneto.github.io/EnvRtype/reference/scan_untested_envs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(maizeG)
met <- sim_met(maizeG, n_env = 10, min.cor = 0.2, max.cor = 0.8,
               min.h2 = 0.3, max.h2 = 0.7, seed = 1)
head(met$data)
env_cor(met, "realised")

## Pair with a simulated envirome that recovers 70% of the same C
W <- sim_W(met$C_env, noise = 0.3, n_var = 50)
} # }
```
