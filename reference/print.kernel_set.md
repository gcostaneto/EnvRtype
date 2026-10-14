# Print a kernel set

Print a kernel set

## Usage

``` r
# S3 method for class 'kernel_set'
print(x, ...)
```

## Arguments

- x:

  object of class `"kernel_set"` from
  [`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md).

- ...:

  ignored.

## Value

Returns `x` invisibly. Called for its side effect of printing a one-line
header and a per-kernel table (kernel name, retained rank, mean diagonal
and rank as a percentage of the number of observations) to the console.
