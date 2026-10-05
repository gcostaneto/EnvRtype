# Scan many sites in memory-bounded chunks

\[scan_untested_envs()\] allocates an \`ndraws x g x m\` array, which is
infeasible for raster-sized \`m\`. This wrapper splits \`W_new\` into
chunks, scans each, and stitches the per-site results. Quantities that
are defined ACROSS new sites – \`vcov\` and \`cor\` – cannot be
recovered from chunks and are returned as \`NULL\` with a warning if
requested.

## Usage

``` r
grid_scan(
  object,
  W_train,
  W_new,
  K_G = NULL,
  chunk_size = 500,
  max_draws = 200,
  verbose = TRUE,
  ...
)
```

## Arguments

- object, W_train, K_G:

  as in \[scan_untested_envs()\].

- W_new:

  covariates for every grid cell / site (m x p).

- chunk_size:

  sites per chunk.

- max_draws:

  MCMC draws retained per chunk. Lower than the scan default because
  mapping needs the posterior mean far more than a tight interval.

- verbose:

  report chunk progress.

- ...:

  passed to \[scan_untested_envs()\].

## Value

list with \`yHat_matrix\`, \`yHat_lower\`, \`yHat_upper\`,
\`yHat_prob_win\`, \`genetic_var\`, \`interpolation\`,
\`exceedance_detail\`, \`meta\`. Class \`scan_untested_envs\` so
downstream functions accept it, but with \`vcov\` and \`cor\` set to
\`NULL\`.
