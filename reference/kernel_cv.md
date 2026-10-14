# k-fold cross-validation predictive ability (out-of-sample)

In-sample yHat overstates accuracy. This masks a fold of phenotypes as
NA, refits, and correlates held-out predictions against observed values.

## Usage

``` r
kernel_cv(
  y,
  data,
  random,
  env,
  gid,
  folds = 5,
  scheme = c("cv1", "cv0"),
  seed = NULL,
  tol = 1e-10,
  scale_kernels = TRUE,
  ...
)
```

## Arguments

- y:

  character. Name of the response column in `data`.

- data:

  data.frame with the response, environment and genotype columns.

- random:

  named list of kernels passed to `kernel_model`.

- env:

  character. Name of the environment column in `data`.

- gid:

  character. Name of the genotype column in `data`.

- folds:

  integer. Number of cross-validation folds. Default 5.

- scheme:

  "cv1" (random cells) or "cv0" (leave-one-environment-out).

- seed:

  integer or NULL. Optional seed for reproducible fold assignment.

- tol:

  numeric. Eigenvalue tolerance passed to the kernel decomposition.

- scale_kernels:

  logical. If `TRUE` (default) kernels are scaled before fitting.

- ...:

  further arguments passed to `kernel_model`.

## Value

A list summarising out-of-sample predictive ability, with components
`scheme` (the validation scheme used), `folds` (number of folds),
`accuracy` (overall correlation between held-out predictions and
observations), `rmse` (overall root-mean-square error), `per_fold` (a
data.frame of per-fold `accuracy` and `rmse`), `fold_id` (integer fold
assignment for each row of `data`), `predicted` (held-out predictions)
and `observed` (the observed response).

## Details

The kernels are decomposed ONCE and reused across all folds: they depend
only on genotypes and environments, never on the phenotype vector, so
masking y does not invalidate them. This is the single largest saving in
the package.
