# Retrieve FAO GAEZ Agro-Ecological Zone Classification for Sites

Extracts the FAO/IIASA Global Agro-Ecological Zone (GAEZ v4) class at a
set of coordinates. By default it reads the GAEZ 57-class raster
**bundled with the package**, so it works offline with no arguments
beyond the coordinates. Companion to
[`get_weather`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md)
and
[`get_soil`](https://gcostaneto.github.io/EnvRtype/reference/get_soil.md).

## Usage

``` r
get_AEZ(
  env.id,
  lat,
  lon,
  file = NULL,
  url = aez_default_url(),
  dir.path = NULL,
  legend = NULL,
  buffer = 0,
  overwrite = FALSE,
  timeout = 1800L,
  verbose = TRUE
)
```

## Arguments

- env.id:

  character vector. Environment/site identifiers.

- lat:

  numeric vector. Latitude in decimal degrees (WGS84), same length as
  `env.id`.

- lon:

  numeric vector. Longitude in decimal degrees (WGS84), same length as
  `env.id`.

- file:

  character. Path to a GAEZ GeoTIFF already on disk. If `NULL`
  (default), `get_AEZ` uses the raster **shipped with the package**
  (`system.file("extdata", "aez_v9v2red_5m_CRUTS32_Hist_8110_100_avg.tif", package = "EnvRtype")`).
  Only if that bundled file is unavailable does it fall back to
  downloading from `url`. Supply `file` to use a different GAEZ release
  or scenario.

- url:

  character. Download URL for a GAEZ raster, used only when `file` is
  `NULL` *and* the bundled raster is missing. Defaults to
  [`aez_default_url()`](https://gcostaneto.github.io/EnvRtype/reference/aez_default_url.md),
  a best-effort template that has *not* been verified against a live FAO
  server – see Details.

- dir.path:

  character. Where to cache a downloaded raster. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html). Use a persistent
  directory to avoid re-downloading.

- legend:

  data.frame or NULL. Lookup table with at least `code` and `label`
  columns. `NULL` (default) uses
  [`aez_legend()`](https://gcostaneto.github.io/EnvRtype/reference/aez_legend.md).
  Pass `NA` to skip labelling and return raw codes only.

- buffer:

  numeric. If greater than 0, the radius in metres of a circular
  neighbourhood summarised around each point (majority class and its
  share), rather than a single-cell lookup. Default 0.

- overwrite:

  boolean. Re-download even if the cached file exists. Default `FALSE`.

- timeout:

  integer. Seconds allowed for the download. Default 1800; GAEZ rasters
  are large and R's 60-second default will abort them.

- verbose:

  boolean. Progress messages. Default `TRUE`.

## Value

A data.frame with one row per environment:

- `env`:

  environment id

- `lat`, `lon`:

  coordinates as supplied

- `aez_code`:

  raster cell value

- `aez_label`, `thermal`, `moisture`, `arable`:

  from the legend, when one is used

- `majority_share`:

  when `buffer > 0`, the proportion of cells in the neighbourhood
  holding the majority class – a purity measure

The raster path used is attached as the `"source"` attribute.

## Details

**Where the raster comes from.** GAEZ v4 has no documented per-point
API: it publishes whole-globe GeoTIFFs through a JavaScript portal.
`get_AEZ` is therefore a raster-extract tool and resolves its raster in
three steps – an explicit `file`, else the layer bundled under
`inst/extdata`, else a download from `url`. The bundled layer is the
GAEZ v4 57-class AEZ historical baseline (`CRUTS32_Hist_8110`, ~5
arc-min), which makes the default call fully offline and reproducible.

**The default URL is unverified.** It is only used when both `file` and
the bundled raster are unavailable. FAO has reorganised GAEZ download
paths between releases and the portal is JavaScript-rendered, so no
stable path could be confirmed. If a download is triggered and fails:

1.  open <https://gaez.fao.org/pages/data-access-download>,

2.  download the AEZ classification layer manually,

3.  pass the file with `file = "path/to/layer.tif"`.

**Coordinates must be WGS84 lon/lat.** They are reprojected to the
raster's CRS automatically, but only if the raster declares one. A
raster with an undefined CRS triggers an error rather than a silent
mis-extraction.

**Interpreting `buffer`.** A single-cell lookup at ~9 km resolution can
be misleading near a zone boundary. With `buffer > 0` the returned
`majority_share` shows how homogeneous the neighbourhood is; a value
near 1 means the site sits well inside one zone, while 0.4 means the
classification is effectively arbitrary at that location.

## References

FAO & IIASA (2021). *Global Agro-Ecological Zones (GAEZ v4)*. Rome, FAO.
<https://gaez.fao.org>

## See also

[`aez_legend`](https://gcostaneto.github.io/EnvRtype/reference/aez_legend.md),
[`get_weather`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md),
[`get_soil`](https://gcostaneto.github.io/EnvRtype/reference/get_soil.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sites <- data.frame(
  env = c("Ames_IA", "Lincoln_NE", "Piracicaba_BR"),
  lat = c(42.03, 40.81, -22.71),
  lon = c(-93.62, -96.68, -47.63))

## Default: uses the raster bundled with the package (offline)
aez <- get_AEZ(sites$env, sites$lat, sites$lon)
aez

## Check zone purity in a 25 km neighbourhood
aez25 <- get_AEZ(sites$env, sites$lat, sites$lon, buffer = 25000)
aez25[, c("env", "aez_label", "majority_share")]

## Raw codes, no labelling
get_AEZ(sites$env, sites$lat, sites$lon, legend = NA)

## A different GAEZ layer you downloaded yourself
get_AEZ(sites$env, sites$lat, sites$lon, file = "GAEZ_AEZ_scenario.tif")
} # }
```
