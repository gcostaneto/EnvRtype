# Collect Soil Features from SoilGrids (250 m)

Queries the ISRIC SoilGrids REST API for soil properties at given
coordinates.

## Usage

``` r
get_soil(
  env.id = NULL,
  lat = NULL,
  lon = NULL,
  variables.names = "clay",
  depths = c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"),
  stat = c("mean", "q_5", "q_50", "q_95", "uncertainty"),
  wide = TRUE,
  sleep = 12,
  tries = 3L,
  timeout = 60,
  verbose = TRUE
)
```

## Arguments

- env.id:

  vector (character). Environment identifiers.

- lat, lon:

  vector (numeric). Coordinates in decimal degrees (WGS84).

- variables.names:

  vector (character). SoilGrids properties. See Details.

- depths:

  vector (character). Depth intervals to request.

- stat:

  character. Which summary to return: `"mean"` (default), `"q_5"`,
  `"q_50"`, `"q_95"` or `"uncertainty"`.

- wide:

  boolean. Return one row per environment with one column per
  property-depth (default), else long format.

- sleep:

  numeric. Minimum seconds between successive API requests. Minimum 12
  (SoilGrids fair-use limit of 5 requests/min). Network time counts
  toward this interval, and duplicate coordinates are fetched only once.

- tries:

  integer. Retry attempts per failed request.

- timeout:

  numeric. Maximum seconds to wait for a single SoilGrids response
  before the request is aborted and retried. Prevents one stalled
  connection from blocking the whole sequential download. Default 60.

- verbose:

  boolean. Print progress messages.

## Value

A data.frame of soil features.

## Details

Accepted properties: `bdod`, `cec`, `cfvo`, `clay`, `nitrogen`, `ocd`,
`ocs`, `phh2o`, `sand`, `silt`, `soc`, `wv0010`, `wv0033`, `wv1500`.

SoilGrids returns values in mapped units (e.g. clay in g/kg x 10). No
unit conversion is applied; the `unit` column records what was served.

## References

Poggio, L. et al. (2021). SoilGrids 2.0: producing soil information for
the globe with quantified spatial uncertainty. *SOIL* 7, 217-240.

## See also

`get_weather`, `get_spatial`

## Author

Brandon Monier, modified by Germano Costa Neto

## Examples

``` r
if (FALSE) { # \dontrun{
get_soil(env.id = "NM", lat = -13.05, lon = -56.05,
         variables.names = c("clay", "nitrogen"))
} # }
```
