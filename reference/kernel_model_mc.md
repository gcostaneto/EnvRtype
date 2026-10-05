# Multi-chain fit with convergence diagnostics (Gelman-Rubin R-hat, ESS)

Runs independent chains with different seeds and reports R-hat and a
crude ESS for the residual and each kernel-scale parameter. R-hat above
1.01 or a tiny ESS means the chain has not mixed: increase
\`iterations\`.

## Usage

``` r
kernel_model_mc(..., n_chains = 3, seed = 1, tol = 1e-10, scale_kernels = TRUE)
```

## Arguments

- ...:

  arguments passed to `kernel_model` (e.g. `y`, `data`, `random`, `env`,
  `gid`).

- n_chains:

  integer. Number of independent chains. Default 3.

- seed:

  integer. Base seed; chain `c` uses `seed + c`. Default 1.

- tol:

  numeric. Eigenvalue tolerance passed to the kernel decomposition.

- scale_kernels:

  logical. If `TRUE` (default) kernels are scaled before fitting.

## Details

The decomposition is shared across chains.
