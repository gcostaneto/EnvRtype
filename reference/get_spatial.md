# Extract Point Estimates from Raster Files

Extracts point values from a raster at the coordinates of each
environment. This function supersedes `extract_GIS()`, which relied on
the retired raster and sp packages.

## Usage

``` r
get_spatial(
  digital.raster = NULL,
  which.raster.number = NULL,
  env.dataframe = NULL,
  env.id = NULL,
  lat = NULL,
  lng = NULL,
  name.feature = NULL,
  merge = TRUE,
  crs = 4326,
  verbose = TRUE
)
```

## Arguments

- digital.raster:

  SpatRaster, or a file path readable by
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).

- which.raster.number:

  integer. Layer(s) to keep when the input is multi-layer.

- env.dataframe:

  data.frame containing the environment id and coordinates.

- env.id:

  character. Name of the environment-id column in `env.dataframe`.

- lat, lng:

  character. Names of the latitude/longitude columns, or numeric vectors
  of coordinates used directly.

- name.feature:

  character. Name(s) for the extracted layer(s).

- merge:

  boolean. Join the extracted values back onto `env.dataframe`.

- crs:

  numeric or character. CRS of the input coordinates. Default 4326
  (WGS84).

- verbose:

  boolean. Print progress messages.

## Value

A data.frame of extracted values, joined to `env.dataframe` when
`merge = TRUE`.

## Details

Coordinates are validated before extraction: latitude must lie in \[-90,
90\] and longitude in \[-180, 180\]. This catches the common lat/lon
transposition, which otherwise yields silent `NA`s or values from the
wrong hemisphere.

## See also

`get_weather`, `get_soil`

## Author

Germano Costa Neto

## Examples

``` r
# \donttest{
if (requireNamespace("terra", quietly = TRUE) &&
    requireNamespace("sf", quietly = TRUE)) {
## A small synthetic elevation raster covering the example site.
elevation <- terra::rast(nrows = 2, ncols = 2, xmin = -57, xmax = -55,
                         ymin = -14, ymax = -12, crs = "EPSG:4326")
terra::values(elevation) <- c(100, 200, 300, 400)
ref <- data.frame(env = "NM", lat = -13.05, lng = -56.05)
get_spatial(env.dataframe = ref, digital.raster = elevation,
            lat = "lat", lng = "lng", env.id = "env",
            name.feature = "Elevation")
}
#> ---------------------------------------------------------------
#> get_spatial -- extracts point values from raster files
#> Last updated at 09/27/2026
#> ---------------------------------------------------------------
#> Based on sf::st_as_sf() and terra::extract()
#> Environmental units ................ 1
#> Raster features .................... 1
#> Done ............................... Wed Oct 14 03:14:55 2026
#> Total time ......................... 0.01s
#>   env    lat    lng Elevation
#> 1  NM -13.05 -56.05       300
# }
```
