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
# \donttest{
phenology_templates()
#>            crop                                  label Tbase Tupper day.type
#> 1         maize                       Maize (Zea mays)  10.0     30    short
#> 2       sorghum              Sorghum (Sorghum bicolor)  10.0     34    short
#> 3  pearl_millet      Pearl millet (Pennisetum glaucum)  10.0     34    short
#> 4     sugarcane             Sugarcane (Saccharum spp.)  12.0     32  neutral
#> 5         wheat              Wheat (Triticum aestivum)   0.0     26     long
#> 6        barley               Barley (Hordeum vulgare)   0.0     26     long
#> 7           oat                     Oat (Avena sativa)   0.0     25     long
#> 8          rice                    Rice (Oryza sativa)  10.0     30    short
#> 9       soybean                  Soybean (Glycine max)  10.0     30    short
#> 10  common_bean       Common bean (Phaseolus vulgaris)  10.0     30    short
#> 11     chickpea             Chickpea (Cicer arietinum)   0.0     30     long
#> 12    groundnut           Groundnut (Arachis hypogaea)  10.0     32  neutral
#> 13       cowpea             Cowpea (Vigna unguiculata)  11.0     32    short
#> 14       canola Canola / oilseed rape (Brassica napus)   5.0     27     long
#> 15    sunflower          Sunflower (Helianthus annuus)   6.0     30  neutral
#> 16       cotton            Cotton (Gossypium hirsutum)  15.6     32  neutral
#> 17       potato             Potato (Solanum tuberosum)   7.0     30     long
#> 18      cassava            Cassava (Manihot esculenta)  13.0     30  neutral
#> 19    safflower       Safflower (Carthamus tinctorius)   5.0     30     long
#> 20      alfalfa    Alfalfa / lucerne (Medicago sativa)   5.0     30     long
#>    p.opt n_stages final_stage final_TT                scale
#> 1   12.5        8          R6     2700        Abendroth V/R
#> 2   12.5        7          PM     2300         Vanderlip GS
#> 3   12.0        7          PM     1700             Maiti GS
#> 4   12.5        4         MAT     6500           phenophase
#> 5   16.0        7          PM     2400      Zadoks-anchored
#> 6   16.0        7          PM     2100      Zadoks-anchored
#> 7   16.0        7          PM     2250      Zadoks-anchored
#> 8   12.5        7          PM     2300           Counce P/R
#> 9   12.8        7          R8     2400    Fehr-Caviness V/R
#> 10  12.5        7          R9     1500             CIAT V/R
#> 11  13.0        6          PM     1900           phenophase
#> 12  12.5        7          R8     2300              Boote R
#> 13  12.5        6          PM     1450           phenophase
#> 14  16.0        7          PM     2100        BBCH-anchored
#> 15  12.5        7          R9     2000 Schneiter-Miller V/R
#> 16  12.5        6          CO     2100     nodes/phenophase
#> 17  14.0        5         MAT     1900           phenophase
#> 18  12.5        5         MAT     6000           phenophase
#> 19  14.0        7          PM     2200        BBCH-anchored
#> 20  14.0        6          LF     1000        Kalu-Fick MSC
show_phenology("safflower")
#> ---------------------------------------------------------------
#> Safflower (Carthamus tinctorius)
#> ---------------------------------------------------------------
#>   Tbase 5 C   Tupper 30 C   long-day (p.opt 14 h)
#>   scale: BBCH-anchored
#> 
#>     EM     <=      95 C d   (+   95)
#>     ROS    <=     450 C d   (+  355)
#>     STE    <=     800 C d   (+  350)
#>     BUD    <=    1100 C d   (+  300)
#>     FL     <=    1400 C d   (+  300)
#>     SF     <=    1800 C d   (+  400)
#>     PM     <=    2200 C d   (+  400)
#> 
#>   source: Demir (2026) Turk. J. Field Crops 31(2) doi:10.17557/tjfc.1998413; Flemmer et al. (2015) Ann. Appl. Biol. 166, 331-339 (BBCH scale)
#>   NOTE:   CARDINALS ARE FIELD-VERIFIED: Tbase 5 C / Tupper 30 C from a 12-year x 5-sowing-date trial; emergence 92.7-96.3 C d was stable across sowing dates. Stages after EM are BBCH-anchored ESTIMATES and are NOT independently verified -- calibrate before relying on them.
#> ---------------------------------------------------------------
# }
```
