# Truncate the GxE kernel to its structural rank.

The GxE kernel has a flat eigenvalue tail that behaves like an identity
matrix and competes with the residual. Dropping it reduces sigma2_GE
inflation (mean abs bias 0.333 -\> 0.130) at essentially no cost in
cross-validated prediction (CV2 0.678 full vs 0.669 truncated).

## Usage

``` r
truncate_gxe_kernel(random, keep_prop = 0.75, pattern = "GE")
```

## Arguments

- random:

  a named list of kernels (the `random` argument of `kernel_model`).

- keep_prop:

  numeric in (0, 1\]. Proportion of the GxE eigenvalue mass to retain.
  Default 0.75.

- pattern:

  character. Regular expression selecting the GxE kernels to truncate.
  Default `"GE"`.

## Details

CAUTION: truncation distorts the G:GxE ratio – rg moved to 0.72 vs a
true 0.595 in testing. Use only when sigma2_GE specifically matters.

Any cached decomposition is dropped, since the kernel changes.
