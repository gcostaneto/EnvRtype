# Default GAEZ Download URL

Returns the URL template
[`get_AEZ`](https://gcostaneto.github.io/EnvRtype/reference/get_AEZ.md)
uses when no local file is given and the bundled raster is unavailable.

## Usage

``` r
aez_default_url(layer = "aez57")
```

## Arguments

- layer:

  character. Which GAEZ layer. Currently `"aez57"`.

## Value

A single character string.

## Details

**This URL is unverified.** It was constructed from the documented GAEZ
v4 portal structure but could not be tested against a live FAO server,
and FAO has reorganised these paths between releases. It is exposed as a
function, rather than hidden inside
[`get_AEZ()`](https://gcostaneto.github.io/EnvRtype/reference/get_AEZ.md),
so it can be inspected and overridden without editing the package:


      get_AEZ(env, lat, lon, url = "https://.../correct/path.tif")

For reproducible analyses, rely on the bundled raster or download the
layer once and pass `file=`.

## See also

[`get_AEZ`](https://gcostaneto.github.io/EnvRtype/reference/get_AEZ.md)

## Examples

``` r
aez_default_url()
#> [1] "https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/LR/aez/57_class/aez_v9v2red_5m_ENSEMBLE_rcp2p6_2020s.tif"
```
