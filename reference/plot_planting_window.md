# Plot genotype performance across candidate planting dates

Curves of predicted performance against sowing date, with credible
ribbons and out-of-envelope dates marked. Dates classified
\`gross_extrapolation\` are desaturated and shaded, and counted in the
subtitle, because position is reported but must not be read as a
reliability guarantee.

## Usage

``` r
plot_planting_window(
  scan,
  dates = NULL,
  date_col = "date",
  genotypes = NULL,
  group = NULL,
  summarise = c("genotype", "group", "overall"),
  ribbon = TRUE,
  ribbon_max = 6,
  mark_envelope = TRUE,
  envelope_level = c("gross", "mild"),
  engine = c("ggplot2", "base"),
  ylab = "Predicted performance",
  title = NULL
)
```

## Arguments

- scan:

  a \`scan_untested_envs\` object, or a table from
  \[planting_window_table()\].

- dates, date_col, genotypes, group:

  passed to \[planting_window_table()\] when \`scan\` is a scan object.

- summarise:

  one of \`"genotype"\` (one curve per genotype), \`"group"\` (mean
  curve per group, requires \`group\`), or \`"overall"\` (single mean
  curve across all selected genotypes).

- ribbon:

  draw credible ribbons. Switched off automatically when more than
  \`ribbon_max\` curves are drawn, since overlapping ribbons are
  unreadable.

- ribbon_max:

  curve count above which ribbons are suppressed.

- mark_envelope:

  shade and desaturate dates outside the envelope.

- envelope_level:

  which positions to mark: \`"gross"\` (default) or \`"mild"\` (also
  marks \`mild_extrapolation\`).

- engine:

  \`"ggplot2"\` (default) or \`"base"\`.

- ylab, title:

  axis label and title.

## Value

a ggplot object (invisibly for \`engine = "base"\`).
