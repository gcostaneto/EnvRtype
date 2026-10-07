# Simulate an Envirome with a Controllable Environmental Diversity

Generates an environmental covariable matrix (environments in rows,
covariables in columns) whose diversity among environments is set by
explicit knobs: the structure (unstructured, mega-environment clusters,
or a continuous gradient), the redundancy among covariables and their
marginal shape. The induced environment correlation and environmental
kernel are attached, so the result feeds
[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md)
(as `W` or `C`) and
[`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md)
directly. It is the environmental analogue of
[`sim_markers`](https://gcostaneto.github.io/EnvRtype/reference/sim_markers.md).

## Usage

``` r
sim_envirome(
  n_env,
  n_var = 20L,
  structure = c("random", "clusters", "gradient"),
  n_cluster = 3L,
  between.cor = 0.2,
  spread = 1,
  redundancy = 0,
  n_blocks = 1L,
  marginal = c("gaussian", "right-skewed", "heavy-tailed", "bounded"),
  env_names = NULL,
  var_names = NULL,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- n_env:

  integer. Number of environments (rows, \>= 2).

- n_var:

  integer. Number of covariables (columns, \>= 2).

- structure:

  character. `"random"` (default), `"clusters"` or `"gradient"`.

- n_cluster:

  integer. Number of mega-environments when `structure = "clusters"`.

- between.cor:

  numeric in \[0, 1). Similarity between clusters (0 = well separated).

- spread:

  numeric. Within-cluster / along-gradient dispersion; larger values
  make environments within a group more diverse.

- redundancy:

  numeric in \[0, 1). Share of covariables collapsed onto drivers
  (collinearity), lowering effective rank.

- n_blocks:

  integer. Number of covariable drivers when `redundancy > 0`.

- marginal:

  character. Marginal of the covariables: `"gaussian"` (default),
  `"right-skewed"`, `"heavy-tailed"` or `"bounded"`.

- env_names, var_names:

  character. Optional names for rows / columns.

- seed:

  integer. RNG seed. The caller's RNG stream is restored on exit.

- verbose:

  logical. Print a diversity summary. Default `TRUE`.

## Value

A matrix of class `"sim_envirome"` (environments x covariables).
Attributes: `"cluster"` (mega-environment of each environment, when
applicable), `"C_env"` (the induced environment correlation), `"K_E"`
(the environmental kernel) and `"diversity"` (`mean_similarity` off the
environment correlation and `eff_dim`, the effective dimension of the
environmental kernel).

## Details

Each environment gets a latent profile that drives its covariables; the
`structure` controls how those profiles are arranged. `"random"` places
environments independently (maximal diversity), `"clusters"` groups them
into `n_cluster` mega-environments (high within-cluster similarity,
between-cluster similarity set by `between.cor`), and `"gradient"`
orders environments along a latent continuum (neighbours similar,
extremes different), as for a latitude or season gradient. `redundancy`
collapses a share of covariables onto `n_blocks` drivers, lowering the
effective rank, and `marginal` reshapes the covariables away from
Gaussian.

## See also

[`sim_W`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md),
[`sim_markers`](https://gcostaneto.github.io/EnvRtype/reference/sim_markers.md),
[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md),
[`env_kernel`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## High diversity: unstructured environments
E_hi <- sim_envirome(20, n_var = 30, structure = "random", seed = 1)

## Mega-environments: three tight clusters
E_cl <- sim_envirome(20, n_var = 30, structure = "clusters",
                     n_cluster = 3, between.cor = 0.1, seed = 1)
attr(E_cl, "diversity")

## A latitudinal / seasonal gradient
E_gr <- sim_envirome(20, n_var = 30, structure = "gradient", seed = 1)

## Feed into sim_met as the envirome (reaction norm) or as a target C
met <- sim_met(markers = sim_markers(100, 500, seed = 1),
               W = as.matrix(E_cl), seed = 1)
KE  <- list(W = attr(E_cl, "K_E"))

## Contrast G x E under low vs high environmental diversity (genetics fixed)
X <- sim_markers(150, 2000, seed = 1)
g_hi <- cor(sim_met(markers = X, C = attr(E_hi, "C_env"), seed = 1)$truth$g)
g_cl <- cor(sim_met(markers = X, C = attr(E_cl, "C_env"), seed = 1)$truth$g)
c(diverse = mean(g_hi[lower.tri(g_hi)]),   # low genetic cor -> more G x E
  clustered = mean(g_cl[lower.tri(g_cl)]))
} # }
```
