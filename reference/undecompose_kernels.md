# Drop cached decompositions

Use after modifying a kernel in place, to force recomputation.

## Usage

``` r
undecompose_kernels(K)
```

## Arguments

- K:

  a kernel list from
  [`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md).

## Value

`K` with decompositions and metadata removed.
