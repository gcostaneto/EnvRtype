# Predict Trained Genotypes at Untested Environments Without Refitting

Transfers a fitted
[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
to environments that were never planted – new locations, new planting
dates, or both – using only their environmental covariates. No
re-fitting and no phenotypes at the new sites are required.

The environmental main effect estimated at the training sites is
*kriged* onto each new site through a linear environmental kernel, and
the whole posterior is carried draw-by-draw so that credible intervals,
win probabilities and an extrapolation diagnostic are all internally
consistent.

**The transfer.** Covariates are standardised on the training moments,
\$\$\tilde{W} = (W - \bar{W}\_{tr}) \\ / \\ s\_{tr},\$\$ and a linear
(product) environmental kernel is formed from the \\p\\ covariates,
\$\$K\_{tt} = \tilde{W}\_{tr}\tilde{W}\_{tr}^{\top}/p, \quad K\_{nt} =
\tilde{W}\_{n}\tilde{W}\_{tr}^{\top}/p, \quad K\_{nn} =
\tilde{W}\_{n}\tilde{W}\_{n}^{\top}/p.\$\$ The new-site environmental
effects are a kriging (Gaussian-process conditional-mean) map of the
training effects, \$\$H = K\_{nt}\\(K\_{tt} + \lambda I)^{-1}, \qquad
\hat{e}\_{n}^{(d)} = H\\e\_{tr}^{(d)},\$\$ applied to every retained
MCMC draw \\d\\. When a GxE kernel is present its interaction effects
are kriged the same way. The ridge \\\lambda\\ stabilises the inverse;
\\H\\ is *not* constrained to be positive, and large negative row
weights are the signal that a site is being reached by extrapolation
(see `weight_negativity` below).

**The prediction.** For draw \\d\\, genotype \\g\\ and new site \\j\\,
\$\$\hat{y}^{(d)}\_{gj} = \mu + u_g^{(d)} + \hat{e}^{(d)}\_{j} +
\widehat{(ge)}^{(d)}\_{gj},\$\$ with \\\mu\\ the training grand mean.
Point predictions are posterior means over draws and the credible
interval uses the empirical quantiles at \\\alpha=(1-\texttt{level})/2\\
and \\1-\alpha\\. The probability that a genotype is the best at a site
is the fraction of draws in which it attains the row maximum, so it
respects the genotype correlation *within* a draw. Per-site genomic
heritability is reported as \\h^2 = \sigma^2_g / (\sigma^2_g +
\sigma^2_e)\\.

## Usage

``` r
scan_untested_envs(
  object,
  W_train,
  W_new,
  K_G = NULL,
  standardize = TRUE,
  lambda = 1e-06,
  level = 0.95,
  max_draws = 500,
  envelope = c("keep", "flag", "mask", "drop"),
  envelope_level = c("gross", "mild")
)
```

## Arguments

- object:

  a fitted model from
  [`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md),
  fitted with `keep_effects = TRUE`. The per-draw effect chains (not
  posterior means) are required; without them the intervals would be
  fabricated.

- W_train:

  training environmental covariate matrix, environments \\\times\\
  covariates (\\e \times p\\). If it has no `rownames` they are taken to
  be the sorted training environments.

- W_new:

  new-environment covariate matrix, sites \\\times\\ covariates (\\m
  \times p\\). Must have the *same columns* as `W_train`. Unnamed rows
  are labelled `NewEnv1..NewEnvm`.

- K_G:

  optional genomic relationship matrix. Used only to inform the GxE
  transfer; when omitted, GxE is kriged from environmental similarity
  alone (a warning is issued).

- standardize:

  logical. Standardise covariates on the training mean and SD before
  building the kernels. Default `TRUE` and strongly advised when
  covariates are on different scales.

- lambda:

  numeric ridge added to \\K\_{tt}\\ before inversion, as a fraction of
  its mean diagonal. Default `1e-6`. Raise it if the kernel is
  near-singular (highly collinear covariates).

- level:

  credible-interval mass. Default `0.95`.

- max_draws:

  integer. Cap on the number of MCMC draws used; draws are thinned by
  even spacing when the chain is longer. Default `500`.

- envelope:

  what to do about sites that fall outside the training covariate
  envelope: `"keep"` (default, return everything unchanged), `"flag"`
  (return everything, warn), `"mask"` (set flagged predictions to `NA`
  but keep the matrix structure so indexing by environment name still
  resolves), or `"drop"` (remove flagged rows/columns and record them in
  `$excluded`).

- envelope_level:

  which sites count as flagged: `"gross"` (default) or `"mild"`. Only
  `"gross"` separates cleanly in calibration; `"mild"` acts on a weak
  signal and is not recommended as a gate.

## Value

An object of class `"scan_untested_envs"`: a list with

- `yHat_matrix`, `yHat_lower`, `yHat_upper`:

  genotype \\\times\\ new-site matrices of the posterior mean and
  credible bounds.

- `yHat`:

  the same in long form, with `ci_width`.

- `yHat_prob_win`:

  probability each genotype is best at each site.

- `genetic_var`:

  per-site genetic variance, its interval and `genomic_h2`.

- `vcov`, `cor`:

  between-new-site genetic covariance and correlation matrices.

- `interpolation`, `divergence`, `diagnostics`:

  the extrapolation measures described above, one row per new site.

- `exceedance_detail`:

  one row per (site \\\times\\ offending covariate) that fell outside
  the training range.

- `excluded`:

  sites masked or dropped by `envelope`, else `NULL`.

- `meta`:

  dimensions, settings, similarity and kriging weights.

## Details

**Extrapolation diagnostics.** The scan reports *where* each new site
sits relative to the training covariate cloud. Several complementary
measures are returned because no single one is sufficient:

- *Cosine similarity* to the nearest training site, \\s\_{ij} =
  K\_{nt,ij}/\sqrt{K\_{nn,ii}\\K\_{tt,jj}}\\, summarised as
  `max_similarity` and `mean_similarity`.

- A *Mahalanobis-type distance* in the leading principal-component
  subspace (enough PCs to reach 95% variance), \\d_i = \sqrt{\sum_k
  z\_{ik}^2}\\, divided by the 95th percentile of the training distances
  to give `mahalanobis_ratio`.

- A *marginal range test* per covariate (`prop_outside`,
  `max_exceedance_sd`). This is marginal: a site inside every
  covariate's range individually can still occupy a *combination* that
  never occurred, which only the distance and weight measures catch.

- *Kriging-weight geometry*: `weight_negativity` \\= \sum_j
  \min(H\_{ij},0)\\ and the effective number of training sites leaned
  on, \\n\_{\mathrm{eff},i} = (\sum_j \|H\_{ij}\|)^2 / \sum_j
  H\_{ij}^2\\.

- A combined `divergence_index` \\= \min\\\big(1, \max(0,\\
  0.5(1-s\_{\max}) + 0.5\\t)\big)\\ with \\t =
  \mathrm{clamp}((\texttt{mahalanobis\\ratio}-0.5)/2,\\0,\\1)\\.

Sites are labelled `"interpolation"`, `"edge"`, `"mild_extrapolation"`
or `"gross_extrapolation"` from these.

**Position is geometry, not a correctness guarantee.** Calibration
against 60 held-out environments found the envelope measures correlate
only \\r \approx -0.37\\ to \\-0.66\\ with realised accuracy. Treat the
position label as a description of where a site lies, never as a verdict
that the prediction is right or wrong. The model-internal
`ci_width_ratio` was the strongest single accuracy correlate
(\\r=-0.66\\); to know your true accuracy, run leave-environments-out
cross-validation on your own data.

## References

Costa-Neto, G., Fritsche-Neto, R. & Crossa, J. (2021) Nonlinear kernels,
dominance, and envirotyping data increase the accuracy of genome-based
prediction in multi-environment trials. *Heredity* 126, 92-106.

Jarquin, D. et al. (2014) A reaction norm model for genomic selection
using high-dimensional genomic and environmental data. *Theoretical and
Applied Genetics* 127, 595-607.

Cressie, N. (1993) *Statistics for Spatial Data*. Wiley. (kriging as a
Gaussian-process conditional mean.)

## See also

[`kernel_model`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md),
[`get_kernel`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md),
[`grid_scan`](https://gcostaneto.github.io/EnvRtype/reference/grid_scan.md),
[`map_scan`](https://gcostaneto.github.io/EnvRtype/reference/map_scan.md),
[`plot_planting_window`](https://gcostaneto.github.io/EnvRtype/reference/plot_planting_window.md),
[`best_planting_date`](https://gcostaneto.github.io/EnvRtype/reference/best_planting_date.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("maizeWTH"); data("maizeYield"); data("maizeG")

## Training covariates (first 100 days after sowing)
vars <- c("FRUE", "PETP", "GDD", "T2M_MAX")
Wtr  <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
                 var.id = vars, statistic = "mean")

## Fit KEEPING the effect chains -- required by the scan
K   <- get_kernel(K_G = list(G = maizeG),
                  K_E = list(W = env_kernel(env.data = Wtr)[[2]]),
                  data = maizeYield, model = "RNMM",
                  env = "env", gid = "gid", y = "value")
fit <- kernel_model(y = "value", data = maizeYield, random = K,
                    env = "env", gid = "gid", keep_effects = TRUE)

## Predict the trained genotypes at new sites (reuse the training W here)
sc <- scan_untested_envs(fit, W_train = Wtr, W_new = Wtr, K_G = maizeG)
sc                                   # print method: per-site diagnostics
head(sc$yHat)                        # predictions + credible intervals
sc$interpolation[, c("env", "position", "mahalanobis_ratio")]
} # }
```
