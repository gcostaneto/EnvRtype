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
if (FALSE) { # \dontrun{
ref <- data.frame(env = "NM", lat = -13.05, lng = -56.05)
get_spatial(env.dataframe = ref, digital.raster = my_elevation_raster,
            lat = "lat", lng = "lng", env.id = "env",
            name.feature = "Elevation")
} # }
```
