# FAO GAEZ Agro-Ecological Zone Legend

Returns the class table used to decode GAEZ agro-ecological zone raster
values into human-readable zone names.

## Usage

``` r
aez_legend(classes = 57)
```

## Arguments

- classes:

  integer. Which legend to return. Currently only `57` (the GAEZ v4
  57-class AEZ layer) is built in.

## Value

A data.frame with `code`, `label`, `thermal`, `moisture` and `arable`
(logical; `FALSE` for water, urban, barren and the constraint classes).

## Details

**Verify before relying on these labels.** FAO has revised AEZ class
definitions between GAEZ releases. The table here follows the GAEZ v4
documentation, but if your raster came from a different release the
codes may not line up.
[`get_AEZ`](https://gcostaneto.github.io/EnvRtype/reference/get_AEZ.md)
will warn when a raster contains codes that are absent from the legend –
that warning is the signal to check.

## See also

[`get_AEZ`](https://gcostaneto.github.io/EnvRtype/reference/get_AEZ.md)

## Examples

``` r
leg <- aez_legend()
head(leg)
#>   code                        label thermal  moisture arable
#> 1    1  Tropics, lowland; semi-arid Tropics semi-arid   TRUE
#> 2    2  Tropics, lowland; sub-humid Tropics sub-humid   TRUE
#> 3    3      Tropics, lowland; humid Tropics     humid   TRUE
#> 4    4  Tropics, lowland; per-humid Tropics per-humid   TRUE
#> 5    5 Tropics, highland; semi-arid Tropics semi-arid   TRUE
#> 6    6 Tropics, highland; sub-humid Tropics sub-humid   TRUE
table(leg$thermal)
#> 
#>                Arctic / very cold   Barren / very sparse vegetation 
#>                                 1                                 1 
#>                            Boreal                            Desert 
#>                                 8                                 4 
#>           Land with poor drainage Land with severe soil constraints 
#>                                 1                                 1 
#>           Land with steep terrain            No data / unclassified 
#>                                 1                                 1 
#>                    Protected area                        Subtropics 
#>                                 1                                16 
#>                         Temperate                           Tropics 
#>                                12                                 8 
#>                  Urban / built-up                      Water bodies 
#>                                 1                                 1 
```
