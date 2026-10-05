# WorldClim Bioclimatic Variables at Trial Sites

Returns the 19 standard bioclimatic indices (BIO1-BIO19) for each site,
from the 1970-2000 WorldClim 2.1 normals.

## Usage

``` r
get_bioclim(
  env.id = NULL,
  lat = NULL,
  lon = NULL,
  vars = 1:19,
  resolution = c("2.5m", "5m", "10m", "30s"),
  dir.path = NULL,
  verbose = TRUE
)
```

## Arguments

- env.id:

  character. Environment identifiers.

- lat, lon:

  numeric. Site coordinates (decimal degrees).

- vars:

  integer or character. Which indices, e.g. `c(1, 12)` or
  `c("BIO1", "BIO12")`. Default all 19.

- resolution:

  character. `"10m"`, `"5m"`, `"2.5m"` or `"30s"`. Default `"2.5m"`.

- dir.path:

  character. Cache directory.

- verbose:

  boolean.

## Value

data.frame with `env`, `LAT`, `LON` and one column per requested index.
Units follow WorldClim: temperature indices in C, precipitation in mm,
BIO3/BIO4/BIO15 dimensionless or percent.

## What these are, and are not

Bioclim variables are LONG-TERM CLIMATE NORMALS, not the weather of your
trial year. They describe the climate a site is drawn from, so they are
appropriate for characterising a TPE, clustering locations, or as static
site covariates. They are NOT a substitute for
[`get_weather()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md) +
[`env_phenology()`](https://gcostaneto.github.io/EnvRtype/reference/env_phenology.md)
when you want to explain what happened in a particular season: a site's
BIO1 is identical in a drought year and a wet year. Mixing normals and
realised weather in one kernel without saying which is which is a common
and serious error.

## References

Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial
resolution climate surfaces for global land areas. *International
Journal of Climatology* 37(12), 4302-4315.
[doi:10.1002/joc.5086](https://doi.org/10.1002/joc.5086)

O'Donnell, M.S. & Ignizio, D.A. (2012) Bioclimatic predictors for
supporting ecological applications in the conterminous United States.
*U.S. Geological Survey Data Series* 691.

## Examples

``` r
if (FALSE) { # \dontrun{
sites <- data.frame(env = c("PIRA", "SETE"),
                    lat = c(-22.7, -19.4), lon = c(-47.6, -44.2))
bc  <- get_bioclim(env.id = sites$env, lat = sites$lat, lon = sites$lon)
bc2 <- get_bioclim(sites$env, sites$lat, sites$lon, vars = c(1, 12))
} # }
```
