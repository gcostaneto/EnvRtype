# Decompose one kernel element (dense or block-diagonal). Scaling is applied BEFORE decomposition so it matches kernel_model().

Decompose one kernel element (dense or block-diagonal). Scaling is
applied BEFORE decomposition so it matches kernel_model().

## Usage

``` r
.kd_decompose_one(
  k,
  ne = NULL,
  tol = 1e-10,
  scale_kernel = TRUE,
  keep_var = 1,
  verbose = FALSE
)
```
