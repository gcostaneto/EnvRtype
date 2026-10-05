# Elevation at Trial Sites from WorldClim 2.1

Returns elevation (m) for each site. Elevation is a strong covariate for
temperature and radiation, so it is often worth carrying alongside the
weather-derived envirotypes – but note it is STATIC: it adds
between-site structure, never within-season dynamics.

## Usage

``` r
get_elevation(
  env.id = NULL,
  lat = NULL,
  lon = NULL,
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

- resolution:

  character. `"10m"`, `"5m"`, `"2.5m"` or `"30s"`. Default `"2.5m"`. 30s
  is ~9.7 GB – use only if you genuinely need ~1 km cells.

- dir.path:

  character. Cache directory. Default is the same cache used by
  [`get_climate_scenario()`](https://gcostaneto.github.io/EnvRtype/reference/get_climate_scenario.md).

- verbose:

  boolean.

## Value

data.frame with `env`, `LAT`, `LON`, `ALT` (m).

## Details

Sites in the sea or outside the land mask return `NA` with a warning
naming them, rather than silently propagating NA into downstream
kernels.

## References

Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial
resolution climate surfaces for global land areas. *International
Journal of Climatology* 37(12), 4302-4315.
[doi:10.1002/joc.5086](https://doi.org/10.1002/joc.5086)

## Examples

``` r
if (FALSE) { # \dontrun{
sites <- data.frame(env = c("PIRA", "SETE"),
                    lat = c(-22.7, -19.4), lon = c(-47.6, -44.2))
alt <- get_elevation(env.id = sites$env, lat = sites$lat, lon = sites$lon)
} # }
```
