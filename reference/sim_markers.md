# Simulate a Marker Matrix with a Controllable Genetic Diversity

Generates a bi-allelic SNP dosage matrix (lines in rows, markers in
columns) whose genetic diversity is set by explicit knobs: the
minor-allele-frequency spectrum, population structure (\\F\_{ST}\\
across subpopulations) and linkage disequilibrium (redundant marker
blocks). The VanRaden genomic relationship matrix is attached, so the
result feeds
[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md)
(as `markers` or via the kinship) and
[`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md)
directly.

## Usage

``` r
sim_markers(
  n_lines,
  n_marker = 1000L,
  maf = c("uniform", "rare", "common"),
  n_subpop = 1L,
  Fst = 0,
  ld = 0,
  ld_blocks = 1L,
  line_names = NULL,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- n_lines:

  integer. Number of lines (rows, \>= 2).

- n_marker:

  integer. Number of markers (columns).

- maf:

  character. Minor-allele-frequency spectrum: `"uniform"` (default),
  `"rare"` (U-shaped, many low-frequency variants) or `"common"`
  (bell-shaped, mostly intermediate).

- n_subpop:

  integer. Number of subpopulations. Lines are assigned in round-robin
  order.

- Fst:

  numeric in \[0, 1). Between-subpopulation differentiation. `0` gives a
  single panmictic population.

- ld:

  numeric in \[0, 1). Fraction of each block's markers forced to equal
  the block lead marker, creating linkage disequilibrium / redundancy.

- ld_blocks:

  integer. Number of linkage blocks when `ld > 0`.

- line_names:

  character. Optional line names (default `G1..Gn`).

- seed:

  integer. RNG seed. The caller's RNG stream is restored on exit.

- verbose:

  logical. Print a diversity summary. Default `TRUE`.

## Value

An integer matrix of class `"sim_markers"` (lines x markers).
Attributes: `"freq"` (per-marker allele frequency), `"subpop"`
(subpopulation of each line), `"G"` (the VanRaden kinship, positive
definite) and `"diversity"` (`He` expected heterozygosity,
`mean_relatedness` and `eff_dim`, the effective dimension of the
kinship).

## Details

Allele frequencies are drawn per marker and genotypes as
\\x\_{ij}\sim\mathrm{Binomial}(2, p\_{j})\\. Population structure
follows the Balding–Nichols model: within subpopulation \\s\\ the
frequency drifts from the ancestral \\p_j\\ as
\\p\_{sj}\sim\mathrm{Beta}\\\big(p_j\tfrac{1-F\_{ST}}{F\_{ST}},
(1-p_j)\tfrac{1-F\_{ST}}{F\_{ST}}\big)\\, so larger \\F\_{ST}\\ means
more between-group differentiation and lower within-group diversity.
Linkage disequilibrium is emulated by copying each block's lead marker
into a fraction `ld` of the other markers in the block, lowering the
effective rank of the kinship.

## See also

[`sim_met`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md),
[`sim_envirome`](https://gcostaneto.github.io/EnvRtype/reference/sim_envirome.md),
[`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## High diversity: many common markers, one panmictic population
X_hi <- sim_markers(150, 2000, maf = "common", seed = 1)

## Low diversity: fewer, rarer markers
X_lo <- sim_markers(150, 300, maf = "rare", seed = 1)

## Structured: four subpopulations with strong drift
X_st <- sim_markers(150, 2000, n_subpop = 4, Fst = 0.15, seed = 1)
attr(X_st, "diversity")

## Feed straight into sim_met (markers mechanism) or via the kinship
met <- sim_met(markers = X_st, n_env = 10, n.qtl = 50, seed = 1)
KG  <- list(G = attr(X_st, "G"))

## Compare phenotype variance under low vs high genetic diversity
E <- sim_envirome(10, 30, structure = "clusters", seed = 1)
v_lo <- var(sim_met(markers = X_lo, C = attr(E, "C_env"), seed = 1)$data$value)
v_hi <- var(sim_met(markers = X_hi, C = attr(E, "C_env"), seed = 1)$data$value)
c(low = v_lo, high = v_hi)
} # }
```
