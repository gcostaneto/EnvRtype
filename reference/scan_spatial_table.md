# Join scan predictions to coordinates

Join scan predictions to coordinates

## Usage

``` r
scan_spatial_table(
  scan,
  coords,
  lon_col = "lon",
  lat_col = "lat",
  date_col = NULL,
  genotypes = NULL
)
```

## Arguments

- scan:

  a \`scan_untested_envs\` (or \[grid_scan()\]) result.

- coords:

  data.frame with \`env\` matching the scan's new environments, plus
  longitude and latitude columns.

- lon_col, lat_col:

  coordinate column names.

- date_col:

  optional column giving a planting date per site, for (C).

- genotypes:

  optional subset.

## Value

long data.frame: \`gid\`, \`env\`, \`lon\`, \`lat\`, \[\`date\`\],
\`yHat\`, \`lower\`, \`upper\`, \`position\`, \`mahalanobis_ratio\`,
\`prop_outside\`, \`weight_negativity\`, \`prob_win\` (posterior
probability the genotype is the best at that site; \`NA\` for scans that
do not carry it).
