#' Clay content (g/kg) from 5 to 15cm depth for Nairobi, Kenya
#'
#' Soil profile for clay content (g/kg) around Nairobi, downloaded from SoilGrids
#' (https://soilgrids.org/). Shipped as a GeoTIFF under \code{inst/extdata} and read
#' with \pkg{terra}.
#'
#'@docType data
#'
#'@name clay_5_15
#'
#'@examples
#'\dontrun{
#' f <- system.file("extdata", "clay_5_15.tif", package = "EnvRtype")
#' clay <- terra::rast(f)
#' terra::plot(clay)
#'}
#' @format A single-layer GeoTIFF raster read with \code{terra::rast()}.
NULL
