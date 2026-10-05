# Envirotype-informed Kernels for Multi-Environment Genomic Prediction

Builds the genomic, envirotypic and genotype-by-environment kernels
required by
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md).
Given a genomic relationship matrix (`p x p` genotypes) and optionally
an environmental relationship matrix (`q x q` environments), it returns
one kernel per model term, already spectrally decomposed and ready to
sample from.

Two computational engines are available and produce numerically
identical results:

- `engine = "factored"` (default):

  Never allocates the `n x n` kernel. Each kernel is represented by its
  eigenvectors `U` (`n x r`) and eigenvalues `s`, obtained from the
  small `p x p` / `q x q` inputs. Memory is O(n\*r) instead of O(n^2),
  which is what makes large multi-environment trials feasible.

- `engine = "dense"`:

  Builds the explicit `n x n` matrix and attaches it as `$Kernel`.
  Slower and quadratic in memory, but needed if you want to inspect or
  post-process the kernel itself.

Set `verbose = TRUE` to trace each step with timings; `verbose = FALSE`
(default) is entirely silent.

## Usage

``` r
get_kernel(
  K_G = NULL,
  K_E = NULL,
  data = NULL,
  model = "RNMM",
  intercept.random = FALSE,
  engine = c("factored", "dense"),
  env = "env",
  gid = "gid",
  y = "value",
  tol = 1e-10,
  keep_var = 1,
  scale_kernels = TRUE,
  verbose = FALSE
)
```

## Arguments

- K_G:

  list of genomic relationship matrices (`p x p`), named. If `NULL`, an
  identity genotype kernel is used.

- K_E:

  list of environmental relationship matrices (`q x q`), named. If
  `NULL`, benchmark genomic-only models are built and `model` is
  restricted to `MM` / `MDs`.

- data:

  data.frame with the environment, genotype and phenotype columns.

- model:

  character. One of `"MM"`, `"MDs"`, `"EMM"`, `"EMDs"`, `"RNMM"`,
  `"RNMDs"`. See *Model structures*.

- intercept.random:

  logical. Append an identity genotype kernel `KG_Gi`, separating
  marker-predictable genetic variance from genotype-specific deviation
  the markers do not capture. Requires replication (a genotype seen in
  more than one environment) to be identifiable; with one observation
  per genotype it is confounded with the residual and a warning is
  issued. Default `FALSE`.

- engine:

  character. `"factored"` (default) or `"dense"`.

- env, gid, y:

  character. Column names in `data`.

- tol:

  numeric. Eigenvalue cutoff, relative to the largest eigenvalue. Must
  match `kernel_model(tol =)`. Default `1e-10`.

- keep_var:

  numeric in (0, 1\]. Retain only the leading eigenvalues carrying this
  proportion of each kernel's trace. Sampler cost scales with total
  rank, so e.g. `0.99` is often a large speedup for negligible change in
  predictions. Default `1` (keep everything).

- scale_kernels:

  logical. Normalise each kernel to `mean(diag(K)) = 1` before
  decomposing. Must match `kernel_model(scale_kernels =)`. Default
  `TRUE`.

- verbose:

  logical. Narrate each step with timings. Default `FALSE`.

## Value

A named list of kernels, class `"kernel_set"`. Each element is a list
with:

- `Type`:

  `"D"` (dense) or `"BD"` (block-diagonal).

- `decomp`:

  list of blocks, each `list(s, U, nr, deltav, type, pos, scaled, tol)`.

- `diag_mean`:

  `mean(diag(K))`, computed in O(n\*r).

- `Kernel`:

  the `n x n` matrix – only when `engine = "dense"`.

Attributes: `"ne"` (observations per environment), `"kd_meta"`
(decomposition settings, so
[`kernel_model()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
can validate the cache), and `"obs"` (index vectors and level names).

## Model structures

|           |                           |                            |
|-----------|---------------------------|----------------------------|
| **model** | **terms**                 | **kernels**                |
| `MM`      | y = mu + G                | KG_G                       |
| `MDs`     | y = mu + G + GE           | KG_G, KGE_GE               |
| `EMM`     | y = mu + E + G            | KE_W, KG_G                 |
| `EMDs`    | y = mu + E + G + GE       | KE_W, KG_G, KGE_GE         |
| `RNMM`    | y = mu + E + G + GxW      | KE_W, KG_G, KGE_GW         |
| `RNMDs`   | y = mu + E + G + GE + GxW | KE_W, KG_G, KGE_GE, KGE_GW |

`GE` is a block-diagonal (within-environment) interaction; `GxW` is the
reaction-norm interaction built from the environmental kernel.

## Why the factored engine is exact

Every kernel here is a congruence transform of a small matrix. With `ig`
mapping observations to genotypes and `ie` to environments, \$\$K_G =
Z_g K_g Z_g' = K_g\[ig, ig\], \quad K_E = K_e\[ie, ie\]\$\$ so the
eigendecomposition of the `n x n` kernel follows from the small one in
O(p^3) rather than O(n^3). Writing \\K = A A'\\ with \\A =
V\[ig,\]\sqrt{D}\\, the nonzero eigenpairs come from the small Gram
matrix \\A'A\\, giving an orthonormal `U` and exact eigenvalues.

The interaction kernel is an elementwise product of two expansions,
which is the expansion of a Kronecker product: with \\j = (ie-1) n_g +
ig\\, \\K\_{GE} = (K_e \otimes K_g)\[j, j\]\\. On a complete balanced
grid its eigenpairs are closed-form (eigenvalues \\d_e \otimes d_g\\,
eigenvectors a Khatri-Rao product) so no large eigendecomposition is
ever required.

## Original Version

Costa-Neto et al (2021)

EnvRtype v.1.2.3, Sep 2026

## Reusing kernels

Kernels depend only on genotypes and environments, never on the
phenotype. Masking `y` for a cross-validation fold does **not**
invalidate them, so build once and reuse across folds, chains and
models.

## References

Costa-Neto G, Galli G, Carvalho HF, Crossa J, Fritsche-Neto R (2021).
EnvRtype: a software to interplay enviromics and quantitative genomics
in agriculture. *G3*, 11(4), jkab040.

## See also

[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md),
[`decompose_kernels`](https://gcostaneto.github.io/EnvRtype/reference/decompose_kernels.md)

## Author

Germano Costa Neto. Refactored version.

## Examples

``` r
if (FALSE) { # \dontrun{
data("maizeYield"); data("maizeG"); data("maizeWTH")

ECs <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
                var.id = c("FRUE", "PETP", "SRAD", "T2M_MAX"),
                statistic = "mean")
KE <- list(W = env_kernel(env.data = ECs)[[2]])
KG <- list(G = maizeG)

## factored engine (default): no n x n matrix is ever allocated
K <- get_kernel(K_G = KG, K_E = KE, data = maizeYield, model = "RNMM",
                env = "env", gid = "gid", y = "value", verbose = TRUE)

fit <- kernel_model(y = "value", data = maizeYield, random = K,
                    env = "env", gid = "gid",
                    iterations = 5000, burnin = 1000, verbose = TRUE)

## trim the eigenvalue tail: usually a large speedup, tiny accuracy cost
K99 <- get_kernel(K_G = KG, K_E = KE, data = maizeYield, model = "RNMM",
                  keep_var = 0.99)

## separate marker-predictable from residual genetic variance
Ki <- get_kernel(K_G = KG, K_E = KE, data = maizeYield, model = "RNMM",
                 intercept.random = TRUE)
} # }
```
