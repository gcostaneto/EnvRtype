# Merge an environmental covariable matrix onto a phenotypic table

Expands an environment-level **W** matrix to the \\n\\ observations of a
phenotypic data.frame, so that each row of the output corresponds to one
genotype-in-environment record.

## Usage

``` r
env_expand(env.data, df.pheno, env.id = "env", skip = 3, verbose = TRUE)
```

## Arguments

- env.data:

  matrix. Environment x covariable matrix with environments as rownames.

- df.pheno:

  data.frame. Phenotypic records containing the environment column.

- env.id:

  character. Name of the environment column in `df.pheno`.

- skip:

  integer. Number of leading columns of the merged table to drop (the
  phenotypic id columns). Default 3.

- verbose:

  boolean. If `TRUE` (default) prints a progress banner.

## Value

A numeric matrix of \\n\\ rows (observations) by \\k\\ covariables.

## See also

`W_matrix`, `env_kernel`

## Examples

``` r
# \donttest{
data("maizeYield"); data("maizeWTH")
W <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
              var.id = c("T2M", "PRECTOT"), statistic = "mean")
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W
Wn <- env_expand(env.data = W, df.pheno = maizeYield, env.id = "env")
#> ---------------------------------------------------------------
#> env_expand -- expands W to phenotypic observations
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - merging covariables onto phenotypic records
dim(Wn)
#> [1] 750   2
# }
```
