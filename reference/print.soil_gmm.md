# Print a soil_gmm object

Compact summary of a fitted soil mixture model: the selected model and
cluster count, the BIC gap to the runner-up, cluster sizes, the number
of uncertain sites, and the per-cluster profile means in original units.

The BIC gap is reported because a point estimate of the zone count is
misleading on its own. Differences are graded weak (\< 2), positive
(2-6), strong (6-10) and very strong (\> 10). A weak gap means several
partitions fit about equally well.

## Usage

``` r
# S3 method for class 'soil_gmm'
print(x, ...)
```

## Arguments

- x:

  an object of class `soil_gmm`.

- ...:

  ignored.

## Value

`x`, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- soil_classification(soil_grid, risk.vars = "soc")
print(fit)
} # }
```
