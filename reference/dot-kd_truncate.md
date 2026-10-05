# Truncate a decomposition to the leading eigenvalues carrying \`keep_var\` of the trace. The Gibbs inner loop is O(T \* k \* n \* r), so dropping the flat eigenvalue tail speeds up EVERY iteration of EVERY fit.

Truncate a decomposition to the leading eigenvalues carrying
\`keep_var\` of the trace. The Gibbs inner loop is O(T \* k \* n \* r),
so dropping the flat eigenvalue tail speeds up EVERY iteration of EVERY
fit.

## Usage

``` r
.kd_truncate(d, keep_var)
```
