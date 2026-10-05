# Easily Building of Environmental Relatedness Kernels

Returns environmental kinships for reaction-norm models. The output is a
list containing `varCov`, the relatedness among environmental
covariables, and `envCov`, the relatedness among environments. Linear
(GBLUP-like), nonlinear Gaussian and arc-cosine deep kernels are
supported.

## Usage

``` r
env_kernel(
  env.data,
  Y = NULL,
  is.scaled = TRUE,
  sd.tol = 10,
  digits = 5,
  tol = 0.001,
  merge = FALSE,
  Z_E = NULL,
  stages = NULL,
  env.id = "env",
  gaussian = FALSE,
  h.gaussian = NULL,
  deep.kernel = FALSE,
  deep.layers = 1,
  QC = FALSE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  matrix. Environmental variables (or markers) per environment (or per
  genotype-environment combination), typically a `W_matrix` output. An
  envirotype frequency matrix from
  [`T_matrix`](https://gcostaneto.github.io/EnvRtype/reference/T_matrix.md)
  is accepted identically (rows = environments); pass `is.scaled = TRUE`
  since `T_matrix` rows are already normalised compositions.

- Y:

  data.frame. Phenotypic data set containing environment id, genotype id
  and trait value. Only used when `merge = TRUE`.

- is.scaled:

  boolean. If the environmental data are already mean-centred and scaled
  (default `TRUE`), assuming \\x \sim N(0,1)\\.

- sd.tol:

  numeric. Maximum standard deviation tolerated in quality control.
  Columns above this value are eliminated. Default 10.

- digits:

  numeric. Number of digits used for rounding. Default 5.

- tol:

  numeric. Numerical tolerance. Default 1E-3.

- merge:

  boolean. If `TRUE`, the environmental covariables are merged with `Y`
  to build an \\n \times n\\ kernel.

- Z_E:

  matrix. Model matrix for environments, used when `merge = TRUE`. If
  `NULL` it is built from `Y`.

- stages:

  vector (character). Names of each stage or time interval. If not
  `NULL`, one kernel is produced per development stage.

- env.id:

  character. Identification of the environment column. Default `'env'`.

- gaussian:

  boolean. If `TRUE`, uses the Gaussian kernel parametrisation \\envCov
  = exp(-h \cdot d / q)\\.

- h.gaussian:

  numeric. Bandwidth \\h\\ used when `gaussian = TRUE`.

- deep.kernel:

  boolean. If `TRUE`, an arc-cosine deep kernel is built from the
  environmental data instead of the linear/Gaussian kernel. Default
  `FALSE`.

- deep.layers:

  integer. Number of hidden layers when `deep.kernel = TRUE` (\>= 1).

- QC:

  boolean. If `TRUE`, applies quality control based on `sd.tol`
  regardless of `is.scaled`. Default `FALSE`.

- verbose:

  boolean. If `TRUE` (default) prints quality-control messages.

## Value

A list with two elements: `varCov`, the covariable x covariable
relatedness, and `envCov`, the environment x environment relatedness.
When `stages` is supplied, each element is itself a named list of
kernels, one per stage.

## Details

Three kernel methods are available.

**Linear (GB).** \\K = WW' / \mathrm{trace}(WW')/n\\, the environmental
analogue of the GBLUP kernel of VanRaden (2008).

**Gaussian (GK).** \\K = \exp(-h\\d/q)\\, with \\d\\ the squared
Euclidean distance and \\q\\ its median, so the bandwidth is scale-free.

**Arc-cosine deep (DK).** Emulates a deep neural network with
`deep.layers` hidden layers (Cuevas et al. 2019); no bandwidth parameter
is required.

Note that quality control is now controlled by the explicit `QC`
argument. In earlier versions it ran only when `is.scaled = FALSE`, so
with the default `is.scaled = TRUE` no filtering happened despite
`sd.tol` being documented as active.

## References

Cuevas J. et al. (2019). Deep kernel for genomic and near infrared
predictions in multi-environment breeding trials. *G3* 9(9), 2913-2924.

Costa-Neto G., Fritsche-Neto R., Crossa J. (2021). Nonlinear kernels,
dominance, and envirotyping data increase the accuracy of genome-based
prediction in multi-environment trials. *Heredity* 126, 92-106.

## See also

`W_matrix`, `get_kernel`, `env_cluster`,
[`T_matrix`](https://gcostaneto.github.io/EnvRtype/reference/T_matrix.md)

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
data('maizeYield'); data("maizeWTH")

## Environmental covariable matrix (environments x covariables)
W.cov <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
                  var.id = c("T2M", "T2M_MAX", "PRECTOT"), statistic = "mean")
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W

## 1. Linear (co)variance kernel
K.lin <- env_kernel(env.data = W.cov, gaussian = FALSE)
#> ---------------------------------------------------------------
#> env_kernel -- builds environmental relatedness kernels
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing a single environmental kernel
round(K.lin$envCov, 2)
#>       NM    SO    PM    IP    SE
#> NM  0.23 -0.11 -0.23  0.32 -0.21
#> SO -0.11  2.68 -0.60 -0.63 -1.34
#> PM -0.23 -0.60  0.46 -0.25  0.62
#> IP  0.32 -0.63 -0.25  0.63 -0.07
#> SE -0.21 -1.34  0.62 -0.07  0.99

## 2. Nonlinear Gaussian kernel
K.gau <- env_kernel(env.data = W.cov, gaussian = TRUE)
#> ---------------------------------------------------------------
#> env_kernel -- builds environmental relatedness kernels
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing a single environmental kernel
round(K.gau$envCov, 2)
#>      NM   SO   PM   IP   SE
#> NM 1.00 0.14 0.49 0.87 0.36
#> SO 0.14 1.00 0.07 0.06 0.02
#> PM 0.49 0.07 1.00 0.37 0.88
#> IP 0.87 0.06 0.37 1.00 0.33
#> SE 0.36 0.02 0.88 0.33 1.00

## 3. Arc-cosine deep kernel with 2 hidden layers
K.deep <- env_kernel(env.data = W.cov, deep.kernel = TRUE, deep.layers = 2)
#> ---------------------------------------------------------------
#> env_kernel -- builds environmental relatedness kernels
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing a single environmental kernel
round(K.deep$envCov, 2)
#>      NM   SO   PM   IP   SE
#> NM 0.55 0.85 0.27 0.80 0.44
#> SO 0.85 6.44 0.99 1.18 1.29
#> PM 0.27 0.99 1.11 0.49 1.52
#> IP 0.80 1.18 0.49 1.52 0.89
#> SE 0.44 1.29 1.52 0.89 2.38

## 4. One kernel per development stage
stages   <- c('VE', 'V1_V6', 'V6_VT', 'VT_R1')
interval <- c(0, 7, 30, 65, 90)
W.stage  <- W_matrix(env.data = maizeWTH, var.id = c('FRUE', 'PETP'),
                     by.interval = TRUE, time.window = interval,
                     names.window = stages)
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W
K.stage <- env_kernel(env.data = W.stage, stages = stages, gaussian = TRUE)
#> ---------------------------------------------------------------
#> env_kernel -- builds environmental relatedness kernels
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing one kernel per developmental stage
names(K.stage$envCov)
#> [1] "VE"    "V1_V6" "V6_VT" "VT_R1"

## 5. Feeding a reaction-norm model
KE <- list(W = env_kernel(env.data = W.cov)$envCov)
#> ---------------------------------------------------------------
#> env_kernel -- builds environmental relatedness kernels
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing a single environmental kernel
# }
```
