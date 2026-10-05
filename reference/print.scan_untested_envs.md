# Print a Scan of Untested Environments

Compact console summary of a
[`scan_untested_envs`](https://gcostaneto.github.io/EnvRtype/reference/scan_untested_envs.md)
result: the per-site prediction and its extrapolation diagnostics, a
tally of position labels, the covariates that fell outside the training
range, and a reminder of how the geometric measures relate to calibrated
accuracy.

## Usage

``` r
# S3 method for class 'scan_untested_envs'
print(x, ...)
```

## Arguments

- x:

  an object of class `"scan_untested_envs"`.

- ...:

  ignored; present for S3 compatibility.

## Value

`x`, invisibly.

## See also

[`scan_untested_envs`](https://gcostaneto.github.io/EnvRtype/reference/scan_untested_envs.md)
