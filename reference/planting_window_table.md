# Tabulate scan predictions against planting date

Reshapes a \`scan_untested_envs()\` result into one row per genotype x
planting date, carrying the credible bounds and the envelope position.
Useful on its own when you want the numbers rather than a plot.

## Usage

``` r
planting_window_table(
  scan,
  dates,
  date_col = "date",
  genotypes = NULL,
  group = NULL
)
```

## Arguments

- scan:

  a \`scan_untested_envs\` object whose new "environments" are candidate
  planting dates.

- dates:

  either (a) a data.frame with a column \`env\` matching the scan's
  new-environment names and a column giving the date, or (b) a vector of
  dates in the same order as the scan's new environments.

- date_col:

  name of the date column when \`dates\` is a data.frame.

- genotypes:

  optional character vector to subset genotypes.

- group:

  optional named vector or data.frame (\`gid\`, \`group\`) assigning
  genotypes to groups, for group-mean curves.

## Value

data.frame with \`gid\`, \`env\`, \`date\`, \`yHat\`, \`lower\`,
\`upper\`, \`ci_width\`, \`position\`, \`mahalanobis_ratio\`,
\`prop_outside\`, \`weight_negativity\`, and \`group\` when supplied.
