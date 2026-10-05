# Delineate Soil Zones with a Gaussian Mixture Model

Clusters soil profiles returned by `get_soil` into a small number of
zones using a Gaussian mixture model, selecting both the number of
clusters and the covariance structure by BIC. Returns hard cluster
labels, the full matrix of posterior membership probabilities, and a
per-site uncertainty measure.

Unlike k-means, a mixture model gives every site a probability of
belonging to every zone. A site on a boundary is genuinely ambiguous,
and that ambiguity is reported rather than hidden behind a hard label.

## Usage

``` r
soil_classification(
  soil.data,
  env.id = "env",
  variables = NULL,
  depths = NULL,
  G = 1:9,
  model.names = NULL,
  scale.data = TRUE,
  clr.texture = TRUE,
  convert.units = TRUE,
  use.pca = FALSE,
  pca.var = 0.95,
  risk.vars = NULL,
  risk.direction = c("high", "low"),
  uncertainty.threshold = 0.7,
  seed = 1,
  verbose = TRUE
)
```

## Arguments

- soil.data:

  data.frame. Wide-format output of `get_soil(wide = TRUE)`, one row per
  site, with columns named `property|depth` (for example
  `"clay|0_5cm"`).

- env.id:

  character. Name of the site-id column. Default `"env"`.

- variables:

  character. Soil properties to cluster on, for example
  `c("clay", "sand", "soc")`. `NULL` (default) uses every property
  present. Filtering is strongly recommended – see **Choosing
  variables** below.

- depths:

  character. Depth layers to use, for example `c("0_5cm", "30_60cm")`.
  `NULL` (default) uses all layers present. Adjacent SoilGrids layers
  are highly correlated, so using all six rarely helps.

- G:

  integer vector. Candidate numbers of clusters. Default `1:9`. `G = 1`
  is included deliberately: if BIC selects it, the data provide no
  evidence for more than one zone, which is a legitimate and useful
  finding.

- model.names:

  character. `mclust` covariance families to try, for example `"EII"`
  (spherical, equal volume) or `"VVV"` (ellipsoidal, all free). `NULL`
  (default) lets `mclust` try every family it can fit. Restrict this
  when sites are few relative to variables.

- scale.data:

  logical. Standardise variables to zero mean and unit variance before
  fitting. Default `TRUE`, and strongly recommended: soil properties
  differ by orders of magnitude, and an unscaled fit is dominated by
  whichever variable happens to have the largest units.

- clr.texture:

  logical. Apply an isometric log-ratio (ilr) transform to texture
  fractions (clay/sand/silt) at each depth. Default `TRUE`.

- convert.units:

  logical. Apply SoilGrids conversion factors, turning mapped integers
  into conventional units. Requires `.SG_FACTORS` from
  `water_balance.R`. Default `TRUE`.

- use.pca:

  logical. Reduce the feature matrix to principal components before
  fitting. Default `FALSE`. Useful when variables are many and
  collinear.

- pca.var:

  numeric in (0, 1\]. Proportion of variance the retained components
  must explain when `use.pca = TRUE`. Default `0.95`.

- risk.vars:

  character. Variables used to order the cluster labels, for example
  `"soc"`. `NULL` (default) orders by the mean of all scaled features.
  See **Label ordering** below.

- risk.direction:

  `"high"` or `"low"`. Whether high values of `risk.vars` should sort
  last (`"high"`, the default, so cluster 1 is the lowest) or first.

- uncertainty.threshold:

  numeric in (0, 1). A site whose maximum posterior probability falls
  below this is flagged `uncertain`. Default `0.7`.

- seed:

  integer or `NULL`. Seed for reproducible EM starts. Default `1`.
  Cluster *labels* are seed-invariant because of the ordering step,
  though the underlying fit may vary.

- verbose:

  logical. Print progress. Default `TRUE`.

## Value

An object of class `soil_gmm`, a list with:

- `classification`:

  data.frame, one row per input site in input order: `env`, `cluster`,
  `max_posterior`, `uncertainty` (= 1 - max_posterior), and logical
  `uncertain`.

- `posterior`:

  matrix of posterior membership probabilities, sites x clusters, each
  row summing to 1.

- `profiles`:

  data.frame of per-cluster means in original units, plus `n_sites`.

- `model`:

  list with `G`, `modelName`, `bic`, `loglik`.

- `BIC`:

  the full BIC table over \\G\\ and covariance family.

- `features`:

  the feature matrix actually clustered.

- `pca`:

  the `prcomp` object, or `NULL`.

- `settings`:

  the call arguments, for reproducibility.

## Details

**What the model does.** A Gaussian mixture assumes the sites are drawn
from \\G\\ multivariate normal components with unknown means,
covariances and mixing weights, fitted by EM. `mclust` additionally
searches over constrained covariance families (spherical, diagonal,
ellipsoidal; equal or varying across components), and BIC selects both
\\G\\ and the family in one criterion. See Fraley & Raftery (2002) and
Scrucca et al. (2016).

**Compositional data.** Texture fractions sum to a constant. Feeding
them to a Gaussian untransformed is a category error: the closure
constraint forces spurious negative correlations. The default
`clr.texture = TRUE` applies an *isometric* log-ratio transform per
depth, mapping \\D\\ parts to \\D-1\\ orthonormal coordinates. A plain
clr transform is deliberately not used for fitting: its coordinates sum
to zero, leaving the covariance singular and causing BIC to select the
wrong number of clusters. Turn this off only if your variables are not
compositional.

**Label ordering.** Mixture component labels are arbitrary and permute
between runs. Clusters are therefore renumbered by `risk.vars` so that
cluster 1 is always the low end and cluster \\G\\ the high end. Without
this, "cluster 3" would mean something different on every run and
results could not be compared or cached.

**Choosing variables.** A mixture estimates a mean vector and a
covariance per component, so the parameter count grows quickly with the
number of features. Passing every property at every depth gives 20+
highly correlated columns and will overwhelm a modest number of sites.
The function warns when `n < 10 * p`. Prefer a handful of properties at
two contrasting depths, or set `use.pca = TRUE`.

**Sample size.** This is the single most common way to get a meaningless
result. Clustering a handful of sites produces a model in which every
site owns its own component: the output reports `max_posterior = 1` and
`uncertainty = 0` for every row, which looks like total confidence and
means nothing. Zone delineation samples a grid of points across the area
of interest – typically hundreds – not one point per farm. Check \\n\\
against \\p\\ before believing any output.

**Reading the BIC gap.** [`print()`](https://rdrr.io/r/base/print.html)
reports the BIC difference between the selected model and the runner-up,
labelled weak (\< 2), positive (2-6), strong (6-10) or very strong (\>
10) following the usual convention for BIC differences. A weak gap means
several partitions fit about equally well and the zone count should not
be presented as settled.

**Relationship to the source paper.** The paper behind this function
(Zhang et al., doi:10.34133/ehs.0402) could not be read: the publisher
blocks automated access. This implements the standard GMM zoning
workflow its title describes – log-ratio transform, standardise,
BIC-selected mixture, posterior-based uncertainty – but the paper's
exact preprocessing, covariance families and thresholds are unverified.
Treat this as a methodologically conventional starting point, not a
reproduction.

**What this cannot do.** SoilGrids contains no cadmium, and no heavy
metal of any kind. The source paper's subject is cadmium risk; this
function can only cluster the properties SoilGrids actually provides
(texture, organic carbon, pH, CEC, water retention). Any "risk"
interpretation rests entirely on whichever proxies you pass to
`risk.vars`, and the paper's machine learning step – predicting the
contaminant before zoning – has no counterpart here. Do not present
output from this function as a contaminant risk map.

## Output granularity

The function classifies *rows*, one output row per input row, in input
order. If each row is a grid point, you get a cluster per grid point; a
farm spanning several rows may well split across zones, and that split
is usually the informative result. Join back to your own table on `env`.

## References

Zhang, F., Jia, Z., Wu, S., Chen, C., Chen, X., Zheng, C., Xu, M. (2025)
Machine Learning and Gaussian Mixture Model for Delineating Soil Cadmium
Risk Zones. *Ecosystem Health and Sustainability*.
[doi:10.34133/ehs.0402](https://doi.org/10.34133/ehs.0402)

Fraley, C. & Raftery, A.E. (2002) Model-based clustering, discriminant
analysis and density estimation. *Journal of the American Statistical
Association* 97(458), 611-631.
[doi:10.1198/016214502760047131](https://doi.org/10.1198/016214502760047131)

Scrucca, L., Fop, M., Murphy, T.B. & Raftery, A.E. (2016) mclust 5:
clustering, classification and density estimation using Gaussian finite
mixture models. *The R Journal* 8(1), 205-233.
[doi:10.32614/RJ-2016-021](https://doi.org/10.32614/RJ-2016-021)

Aitchison, J. (1982) The statistical analysis of compositional data.
*Journal of the Royal Statistical Society B* 44(2), 139-177.

Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. & Barcelo-Vidal,
C. (2003) Isometric logratio transformations for compositional data
analysis. *Mathematical Geology* 35(3), 279-300.
[doi:10.1023/A:1023818214614](https://doi.org/10.1023/A%3A1023818214614)

Schwarz, G. (1978) Estimating the dimension of a model. *Annals of
Statistics* 6(2), 461-464.
[doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)

Poggio, L. et al. (2021) SoilGrids 2.0. *SOIL* 7, 217-240.
[doi:10.5194/soil-7-217-2021](https://doi.org/10.5194/soil-7-217-2021)

## See also

`get_soil` for the input, `water_balance` for per-site indices that need
no cross-site sample, `Mclust` for the underlying model.

## Examples

``` r
if (FALSE) { # \dontrun{
library(mclust)

## ---------------------------------------------------------------
## 1. Minimal use
## ---------------------------------------------------------------
soil <- get_soil(env.id = c("NM", "SO", "IP"),
                 lat = c(-13.05, -12.32, -21.98),
                 lon = c(-56.08, -55.71, -47.88),
                 variables.names = c("clay", "sand", "silt", "soc", "phh2o"))

fit <- soil_classification(soil)
fit                                  # summary, BIC gap, cluster profiles
head(fit$classification)             # per-site cluster + uncertainty


## ---------------------------------------------------------------
## 2. Realistic zoning: sample a GRID, not one point per farm
## ---------------------------------------------------------------
## A mixture cannot be identified from a handful of rows. Sample the
## area of interest, here 40 points around each of 5 locations.
sites <- data.frame(
  env = c("SOR", "LON", "SLG", "PAL", "BAR"),
  lat = c(-12.5453, -23.3045, -28.4083, -10.1689, -12.1530),
  lon = c(-55.7113, -51.1696, -54.9608, -48.3317, -44.9900))

g <- expand.grid(k = 1:40, i = seq_len(nrow(sites)))
pts <- data.frame(
  env  = paste0(sites$env[g$i], "_", sprintf("%03d", g$k)),
  farm = sites$env[g$i],
  lat  = sites$lat[g$i] + runif(nrow(g), -0.05, 0.05),
  lon  = sites$lon[g$i] + runif(nrow(g), -0.05, 0.05))

soil_grid <- get_soil(env.id = pts$env, lat = pts$lat, lon = pts$lon,
                      variables.names = c("clay", "sand", "silt",
                                          "soc", "phh2o"))

fit <- soil_classification(soil_grid,
                           variables = c("clay", "sand", "silt",
                                         "soc", "phh2o"),
                           depths    = c("0_5cm", "30_60cm"),
                           risk.vars = "soc")
fit

## Do farms split across zones? Usually the interesting question.
table(farm = pts$farm, cluster = fit$classification$cluster)


## ---------------------------------------------------------------
## 3. Posterior probabilities and ambiguous sites
## ---------------------------------------------------------------
round(head(fit$posterior), 3)        # membership in EVERY zone

## Sites the model is least sure about -- these sit on a boundary.
amb <- fit$classification[fit$classification$uncertain, ]
amb[order(amb$max_posterior), ]

## Tighten or relax the flag without refitting the model:
fit2 <- soil_classification(soil_grid, uncertainty.threshold = 0.9)
sum(fit2$classification$uncertain)


## ---------------------------------------------------------------
## 4. Ordering zones so cluster 1 always means the same thing
## ---------------------------------------------------------------
lo  <- soil_classification(soil_grid, risk.vars = "soc",
                           risk.direction = "high")   # 1 = lowest SOC
hi  <- soil_classification(soil_grid, risk.vars = "soc",
                           risk.direction = "low")    # 1 = highest SOC
lo$profiles[, c("cluster", "n_sites")]


## ---------------------------------------------------------------
## 5. Many correlated variables: reduce first
## ---------------------------------------------------------------
fit_pca <- soil_classification(soil_grid, use.pca = TRUE, pca.var = 0.95)
fit_pca$pca                          # how many components were kept

## Or restrict the covariance family when sites are few:
fit_par <- soil_classification(soil_grid, model.names = c("EII", "VII",
                                                          "EEI", "VVI"))


## ---------------------------------------------------------------
## 6. Joining results back to your own table
## ---------------------------------------------------------------
out <- merge(pts, fit$classification, by = "env")
head(out[, c("env", "farm", "lat", "lon", "cluster",
             "max_posterior", "uncertain")])

## Dominant zone per farm, and whether the farm is homogeneous
with(out, table(farm, cluster))


## ---------------------------------------------------------------
## 7. Is the cluster count actually supported?
## ---------------------------------------------------------------
fit$BIC                              # full table over G and family
fit$model$G                          # selected number of zones
## print(fit) reports the gap to the runner-up; a weak gap (< 2) means
## the zone count is not settled, whatever the point estimate says.
} # }
```
