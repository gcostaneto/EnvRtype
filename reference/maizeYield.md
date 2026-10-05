# Maize Phenotypic Data Set

This data set was included in Souza et al. (2017) and Cuevas et al
(2019) is from the Helix Seeds Company (HEL). It consists of grain yied
from 452 maize hybrids obtained by crossing 111 pure lines (inbreds);
the hybrids were evaluated in 2015 at five Brazilian sites (E1-E5). The
experimental design used in each site was a randomized block with two
replicates per hybrid. However, to facilitate the demonstration of
functions, only 150 hybrids per environment are being considered, thus
counting 750 genotype x environment observations. Grain yield data are
mean-centered and scaled.

## Usage

``` r
data(maizeYield)
```

## Format

A data frame with 750 rows and 3 variables:

- env:

  environmental id (factor)

- gid:

  genotypic id (factor)

- value:

  grain yield value of genotype plus genotype by location interaction
  effects in kg ha-1 (numeric)

## Examples

``` r
data(maizeYield); head(maizeYield)
#>   env  gid      value
#> 1  NM G001 -0.6282005
#> 2  NM G002 -2.3118840
#> 3  NM G003 -1.0597583
#> 4  NM G004 -0.6221222
#> 5  NM G005  0.1315846
#> 6  NM G006 -0.6646701
```
