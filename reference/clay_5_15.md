# Clay content (g/kg) from 5 to 15cm depth for Nairobi, Kenya

Soil profile for clay content (g/kg) around Nairobi, downloaded from
SoilGrids (https://soilgrids.org/). Shipped as a GeoTIFF under
`inst/extdata` and read with terra.

## Format

A single-layer GeoTIFF raster read with
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).

## Examples

``` r
if (FALSE) { # \dontrun{
f <- system.file("extdata", "clay_5_15.tif", package = "EnvRtype")
clay <- terra::rast(f)
terra::plot(clay)
} # }
```
