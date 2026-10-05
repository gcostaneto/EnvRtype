# Validated variance components, heritabilities and genetic correlation.

var_total = var(y) (PHENOTYPIC; var(yHat) carries no residual and yields
a NEGATIVE residual by subtraction). var_resid = var_total - (varG +
varE + varGE) h2_entry = varG / (varG + varE/q + varGE/(p\*q) +
var_resid) varG is retained in the denominator; omitting it produced h2
\> 1 even with true inputs (observed 1.025 and 0.983).

## Usage

``` r
varcomp_summary(
  fit,
  method = c("postmean", "gibbs"),
  y = NULL,
  p = NULL,
  q = NULL,
  digits = 4
)
```

## Arguments

- fit:

  a fitted `kernel_model` object (or a compatible list carrying variance
  components).

- method:

  character. `"postmean"` (default) uses posterior-mean effects;
  `"gibbs"` uses the stored Gibbs variance components.

- y:

  numeric. Observed phenotypes used as the total-variance reference.
  Taken from `fit$meta$y` when `NULL`.

- p:

  integer. Number of genotypes. Inferred from `fit$meta$gid` when
  `NULL`.

- q:

  integer. Number of environments. Inferred from `fit$meta$env` when
  `NULL`.

- digits:

  integer. Number of decimal places in the summary table. Default 4.
