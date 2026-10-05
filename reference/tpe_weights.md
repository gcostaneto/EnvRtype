# Target Population of Environments (TPE) Weights

Converts an
[`env_risk_profile()`](https://gcostaneto.github.io/EnvRtype/reference/env_risk_profile.md)
into per-environment-type weights summing to one, so that selection
indices, weighted GxE models and testing network decisions can be made
proportional to how often each environment type actually occurs – rather
than how often it happened to be trialled.

## Usage

``` r
tpe_weights(
  risk,
  stage = NULL,
  var = NULL,
  by = c("class", "site"),
  measure = c("freq", "days"),
  drop.benign = FALSE,
  verbose = TRUE
)
```

## Arguments

- risk:

  An `"env_risk"` object from
  [`env_risk_profile()`](https://gcostaneto.github.io/EnvRtype/reference/env_risk_profile.md).

- stage:

  character. Restrict to one stage (e.g. `"R1"`). If `NULL`, weights are
  computed across all stages jointly.

- var:

  character. Restrict to one variable. If `NULL`, classes from all
  variables are pooled, which is only sensible if they are disjoint.

- by:

  character. `"class"` (default) weights environment TYPES; `"site"`
  weights SITES by how representative they are of the TPE.

- measure:

  character. `"freq"` (probability of occurrence, default) or `"days"`
  (expected days, which weights by intensity).

- drop.benign:

  boolean. Exclude classes named `safe`, `optimal`, `none`,
  `non_limiting`. Default `FALSE`.

- verbose:

  boolean. Print the weight table.

## Value

A data.frame with the grouping column, the raw measure, and `weight`
summing to 1. Attribute `"measure"` records which measure was used.

## Details

**What "site weight" means.** With `by = "site"` the weight is the share
of the TPE's total stress-occurrence mass that a site contributes. A
high weight means the site frequently expresses the stresses of
interest, NOT that it is a good discriminator of genotypes – those are
different criteria, and a site can be highly representative yet have
poor heritability. Combine these weights with per-site repeatability
before cutting a testing location.

**Weights are as good as the record.** They inherit every sampling
limitation of the underlying profile, and they assume the historical
climatology still holds.

## See also

[`env_risk_profile`](https://gcostaneto.github.io/EnvRtype/reference/env_risk_profile.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rp <- env_risk_profile(wth, planting = "10-15", crop = "maize")

## How much of the TPE is each drought class at flowering?
tpe_weights(rp, stage = "R1", var = "PETP")

## Which sites carry the most of the TPE's stress mass?
tpe_weights(rp, var = "PETP", by = "site", drop.benign = TRUE)
} # }
```
