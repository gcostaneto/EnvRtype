# List Built-in Phenology Templates

`phenology_templates()` returns the built-in crop templates as a
data.frame; `show_phenology()` prints them, optionally for one crop with
its stage table, provenance and caveats.

Every template is a published tabulation for a REFERENCE cultivar. They
give a consistent relative scale across sites – which is what fixes the
calendar-day misalignment in
[`W_matrix()`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md)
– but they are not calibrated predictions for your germplasm. Maturity
group alone shifts maize R6 by 300-500 C d. Validate against observed
flowering dates before treating stage labels as phenology.

## Usage

``` r
phenology_templates()

show_phenology(crop = NULL)
```

## Arguments

- crop:

  character. Optional single crop; if omitted, all are listed.

## Value

`phenology_templates()` a data.frame, one row per crop.
`show_phenology()` returns its input invisibly and is called for the
side effect of printing.

## Examples

``` r
if (FALSE) { # \dontrun{
phenology_templates()
show_phenology("safflower")
} # }
```
