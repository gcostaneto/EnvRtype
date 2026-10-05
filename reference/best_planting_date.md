# Best planting date per genotype

Best planting date per genotype

## Usage

``` r
best_planting_date(
  scan,
  dates = NULL,
  date_col = "date",
  genotypes = NULL,
  group = NULL,
  exclude_flagged = TRUE,
  envelope_level = c("gross", "mild")
)
```

## Arguments

- scan, dates, date_col, genotypes, group:

  as in \[planting_window_table()\].

- exclude_flagged:

  drop dates outside the envelope before ranking. Default TRUE: an
  optimum found in a region the training data does not cover is not an
  optimum you can act on.

- envelope_level:

  \`"gross"\` (default) or \`"mild"\`.

## Value

data.frame: \`gid\`, \`best_date\`, \`yHat\`, \`lower\`, \`upper\`,
\`position\`, \`n_dates_considered\`.

## Examples

``` r
## A scan is normally produced by scan_untested_envs(); here we build a
## minimal one by hand -- 3 genotypes across 6 candidate sowing windows --
## so the example is self-contained.
gids <- paste0("G", 1:3)
envs <- paste0("sow", 1:6)
mu <- matrix(c(6, 5.5, 5), 3, 6) +
      matrix(rep(c(0, .4, .7, .6, .2, -.3), each = 3), 3, 6)  # peaks mid-window
dimnames(mu) <- list(gids, envs)
scan <- list(
  yHat_matrix = mu, yHat_lower = mu - 0.5, yHat_upper = mu + 0.5,
  interpolation = data.frame(
    env               = envs,
    position          = c("interpolation", "interpolation", "edge",
                          "interpolation", "mild_extrapolation",
                          "gross_extrapolation"),
    mahalanobis_ratio = c(0.6, 0.8, 1.4, 0.9, 2.1, 3.6),
    prop_outside      = c(0, 0, 0.1, 0, 0.3, 0.7),
    weight_negativity = 0, stringsAsFactors = FALSE))
class(scan) <- "scan_untested_envs"

## Map each scanned environment to a calendar sowing date
dates <- data.frame(env = envs,
                    date = as.Date("2024-04-01") + c(0, 14, 28, 42, 56, 70))

## Best *supported* sowing date per genotype (gross extrapolations dropped)
best_planting_date(scan, dates)
#>   gid  best_date yHat lower upper position n_dates_considered
#> 1  G1 2024-04-29  6.7   6.2   7.2     edge                  5
#> 2  G2 2024-04-29  6.2   5.7   6.7     edge                  5
#> 3  G3 2024-04-29  5.7   5.2   6.2     edge                  5

## Rank every date, including the out-of-envelope ones
best_planting_date(scan, dates, exclude_flagged = FALSE)
#>   gid  best_date yHat lower upper position n_dates_considered
#> 1  G1 2024-04-29  6.7   6.2   7.2     edge                  6
#> 2  G2 2024-04-29  6.2   5.7   6.7     edge                  6
#> 3  G3 2024-04-29  5.7   5.2   6.2     edge                  6

## Stricter: also drop mild extrapolations before ranking
best_planting_date(scan, dates, envelope_level = "mild")
#>   gid  best_date yHat lower upper position n_dates_considered
#> 1  G1 2024-04-29  6.7   6.2   7.2     edge                  4
#> 2  G2 2024-04-29  6.2   5.7   6.7     edge                  4
#> 3  G3 2024-04-29  5.7   5.2   6.2     edge                  4
```
