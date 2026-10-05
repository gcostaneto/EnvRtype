# Copula-based Dependence Indices for an Environmental Covariable Matrix

Performs a copula analysis of an environment x covariable matrix (`W`)
and returns a matrix of *copula indices* that can replace the raw
covariables. Sklar's theorem lets the joint distribution of the
environmental variables be split into (i) the marginal behaviour of each
variable and (ii) the dependence structure linking them (the copula). By
moving to the copula scale, each environment is described by where it
sits in the *multivariate* environmental distribution rather than by raw
physical units.

## Usage

``` r
env_copula(
  W,
  index = c("pobs", "joint", "survival", "kendall", "all"),
  groups = NULL,
  ties.method = "average",
  normalize = TRUE,
  verbose = TRUE
)
```

## Arguments

- W:

  matrix or data.frame. Environmental covariables with environments in
  rows and variables in columns (typically a `W_matrix` output).

- index:

  character. Which copula index to return:

  `"pobs"`

  :   Pseudo-observations \\u\_{ij} = r\_{ij}/(q+1)\\, the marginal
      probability-integral transform. Same dimension as `W`; one column
      per variable.

  `"joint"`

  :   Multivariate empirical copula \\C_q(u_i)\\, the joint
      non-exceedance probability of environment \\i\\: the probability
      that all variables are simultaneously at or below their observed
      levels. One column per block.

  `"survival"`

  :   Survival (joint exceedance) copula: the probability that all
      variables are simultaneously at or above their observed levels.
      Highlights compound-stress environments. One column per block.

  `"kendall"`

  :   Kendall function \\K_C(z) = P(C(U) \le z)\\ evaluated at each
      environment's joint probability. A scalar multivariate return
      level ordering environments from jointly-benign to
      jointly-extreme. One column per block.

  `"all"`

  :   Column-binds `pobs`, `joint`, `survival`, `kendall`.

- groups:

  list or NULL. Optional named list splitting the columns of `W` into
  variable blocks, e.g.
  `list(heat = c("T2M_mean","T2M_MAX_mean"), water = "PRECTOT_mean")`.
  Elements may be column names or column indices. A separate
  multivariate copula is fitted per block, so joint/survival/Kendall
  describe dependence *within* a block. If `NULL` (default) all columns
  form a single block named `"all"`.

- ties.method:

  character. How ranks are computed when there are ties, passed to
  `rank`. Default `"average"`.

- normalize:

  boolean. If `TRUE` (default) the Kendall index is rescaled to
  \\\[0,1\]\\ by its empirical range when non-degenerate; `FALSE` keeps
  raw \\K_C\\.

- verbose:

  boolean. If `TRUE` (default) prints progress messages.

## Value

A numeric matrix with \\q\\ environments in rows (rownames preserved
from `W`) and the requested copula indices in columns. Attributes
`"copula.index"` and `"copula.groups"` record the index produced and the
block structure used.

## Details

**Why rank-based.** All indices are computed from ranks, so they are
invariant to any monotone transformation of the marginals: storing
precipitation in mm or log(mm) gives the same result, and no
distributional assumption is imposed. This is the main practical gain
over a z-scored **W**, where one extreme site can dominate the Euclidean
geometry.

**Pseudo-observations.** For \\q\\ environments and variable \\j\\,
\$\$u\_{ij} = r\_{ij}/(q+1)\$\$ with \\r\_{ij}\\ the rank of environment
\\i\\. Dividing by \\q+1\\ keeps values strictly inside \\(0,1)\\,
avoiding boundary problems in downstream density evaluation.

**Empirical copula.** The joint index is the empirical copula evaluated
at each environment's own pseudo-observation, a nonparametric estimator
of \\C\\ that requires no copula family to be chosen.

**Kendall function.** \\K_C\\ is the CDF of \\Z = C(U)\\, estimated by
the Genest-Rivest empirical construction. It collapses a
\\d\\-dimensional dependence structure into one interpretable number per
environment and underpins multivariate return periods in hydrology.

**Sample-size caveat.** These are rank statistics *across environments*,
so their resolution is bounded by \\q\\. With five environments each
pseudo-observation can only take values in 1/6,...,5/6 and the joint
copula is very coarse; the indices become informative with dozens of
environments. A warning is issued when \\q \< 10\\.

**Scope.** The dependence modelled is that *among environmental
variables across the environment panel*. Geographic coordinates are not
used and no spatial autocorrelation (variogram / distance decay) is
fitted. Use `groups` to control which variables share a copula.

## References

Sklar, A. (1959). Fonctions de repartition a n dimensions et leurs
marges. *Publ. Inst. Statist. Univ. Paris* 8, 229-231.

Genest, C., Rivest, L.-P. (1993). Statistical inference procedures for
bivariate Archimedean copulas. *JASA* 88, 1034-1043.

Salvadori, G., De Michele, C., Kottegoda, N.T., Rosso, R. (2007).
*Extremes in Nature: An Approach Using Copulas*. Springer.

## See also

`W_matrix`, `env_kernel`, `env_cluster`

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
data("maizeWTH")
env.data <- maizeWTH[maizeWTH$daysFromStart < 100, ]

W <- W_matrix(env.data = env.data, var.id = c("T2M", "T2M_MAX", "PRECTOT"),
              statistic = "mean")
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - centring, scaling and quality-controlling W

## 1. Marginal copula scale (pseudo-observations)
round(env_copula(W, index = "pobs"), 3)
#> ---------------------------------------------------------------
#> env_copula -- reparameterises covariables into copula indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: Only 5 environments: copula indices are rank statistics across environments and will be very coarse. Interpret with care.
#> env_copula: 5 environments, 3 covariables, 1 block(s); index = 'pobs'.
#>    cop_u_PRECTOT_mean cop_u_T2M_MAX_mean cop_u_T2M_mean
#> NM              0.667              0.500          0.667
#> SO              0.167              0.833          0.833
#> PM              0.333              0.333          0.333
#> IP              0.833              0.667          0.500
#> SE              0.500              0.167          0.167
#> attr(,"copula.index")
#> [1] "pobs"
#> attr(,"copula.groups")
#> attr(,"copula.groups")$all
#> [1] "PRECTOT_mean" "T2M_MAX_mean" "T2M_mean"    
#> 

## 2. Joint non-exceedance probability per environment
env_copula(W, index = "joint")
#> ---------------------------------------------------------------
#> env_copula -- reparameterises covariables into copula indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: Only 5 environments: copula indices are rank statistics across environments and will be very coarse. Interpret with care.
#> env_copula: 5 environments, 3 covariables, 1 block(s); index = 'joint'.
#>    cop_C_all
#> NM 0.5000000
#> SO 0.1666667
#> PM 0.1666667
#> IP 0.5000000
#> SE 0.1666667
#> attr(,"copula.index")
#> [1] "joint"
#> attr(,"copula.groups")
#> attr(,"copula.groups")$all
#> [1] "PRECTOT_mean" "T2M_MAX_mean" "T2M_mean"    
#> 

## 3. Compound-stress (joint exceedance) index
env_copula(W, index = "survival")
#> ---------------------------------------------------------------
#> env_copula -- reparameterises covariables into copula indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: Only 5 environments: copula indices are rank statistics across environments and will be very coarse. Interpret with care.
#> env_copula: 5 environments, 3 covariables, 1 block(s); index = 'survival'.
#>    cop_surv_all
#> NM    0.1666667
#> SO    0.1666667
#> PM    0.5000000
#> IP    0.1666667
#> SE    0.5000000
#> attr(,"copula.index")
#> [1] "survival"
#> attr(,"copula.groups")
#> attr(,"copula.groups")$all
#> [1] "PRECTOT_mean" "T2M_MAX_mean" "T2M_mean"    
#> 

## 4. Kendall multivariate return level
env_copula(W, index = "kendall")
#> ---------------------------------------------------------------
#> env_copula -- reparameterises covariables into copula indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: Only 5 environments: copula indices are rank statistics across environments and will be very coarse. Interpret with care.
#> env_copula: 5 environments, 3 covariables, 1 block(s); index = 'kendall'.
#>    cop_K_all
#> NM         1
#> SO         0
#> PM         0
#> IP         1
#> SE         0
#> attr(,"copula.index")
#> [1] "kendall"
#> attr(,"copula.groups")
#> attr(,"copula.groups")$all
#> [1] "PRECTOT_mean" "T2M_MAX_mean" "T2M_mean"    
#> 

## 5. Blocks: heat and water handled as separate copulas
env_copula(W, index = "all",
           groups = list(heat  = c("T2M_mean", "T2M_MAX_mean"),
                         water = "PRECTOT_mean"))
#> ---------------------------------------------------------------
#> env_copula -- reparameterises covariables into copula indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: Only 5 environments: copula indices are rank statistics across environments and will be very coarse. Interpret with care.
#> env_copula: 5 environments, 3 covariables, 2 block(s); index = 'all'.
#>    cop_u_PRECTOT_mean cop_u_T2M_MAX_mean cop_u_T2M_mean cop_C_heat cop_C_water
#> NM          0.6666667          0.5000000      0.6666667  0.5000000   0.6666667
#> SO          0.1666667          0.8333333      0.8333333  0.8333333   0.1666667
#> PM          0.3333333          0.3333333      0.3333333  0.3333333   0.3333333
#> IP          0.8333333          0.6666667      0.5000000  0.5000000   0.8333333
#> SE          0.5000000          0.1666667      0.1666667  0.1666667   0.5000000
#>    cop_surv_heat cop_surv_water cop_K_heat cop_K_water
#> NM     0.3333333      0.3333333       0.75        0.75
#> SO     0.1666667      0.8333333       1.00        0.00
#> PM     0.6666667      0.6666667       0.25        0.25
#> IP     0.3333333      0.1666667       0.75        1.00
#> SE     0.8333333      0.5000000       0.00        0.50
#> attr(,"copula.index")
#> [1] "all"
#> attr(,"copula.groups")
#> attr(,"copula.groups")$heat
#> [1] "T2M_mean"     "T2M_MAX_mean"
#> 
#> attr(,"copula.groups")$water
#> [1] "PRECTOT_mean"
#> 

## 6. Straight from W_matrix() via the 'copula' argument
Wc <- W_matrix(env.data = env.data, var.id = c("T2M", "T2M_MAX", "PRECTOT"),
               statistic = "mean", copula = "kendall",
               center = FALSE, scale = FALSE)
#> ---------------------------------------------------------------
#> W_matrix -- builds the environmental covariable (W) matrix
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - summarising weather data into environmental covariables
#>   - reparameterising covariables with copula index
#> ---------------------------------------------------------------
#> env_copula -- reparameterises covariables into copula indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: Only 5 environments: copula indices are rank statistics across environments and will be very coarse. Interpret with care.
#> env_copula: 5 environments, 3 covariables, 1 block(s); index = 'kendall'.
#> W replaced by copula index 'kendall' (1 columns).
#>   - centring, scaling and quality-controlling W

## 7. Copula-scale kernel for reaction-norm models
Kc <- env_kernel(env.data = env_copula(W, index = "pobs"), gaussian = TRUE)
#> ---------------------------------------------------------------
#> env_kernel -- builds environmental relatedness kernels
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#>   - computing a single environmental kernel
#> ---------------------------------------------------------------
#> env_copula -- reparameterises covariables into copula indices
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Warning: Only 5 environments: copula indices are rank statistics across environments and will be very coarse. Interpret with care.
#> env_copula: 5 environments, 3 covariables, 1 block(s); index = 'pobs'.
round(Kc$envCov, 2)
#>      NM   SO   PM   IP   SE
#> NM 1.00 0.37 0.53 0.81 0.37
#> SO 0.37 1.00 0.26 0.22 0.08
#> PM 0.53 0.26 1.00 0.37 0.81
#> IP 0.81 0.22 0.37 1.00 0.30
#> SE 0.37 0.08 0.81 0.30 1.00
# }
```
