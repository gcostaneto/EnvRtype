# Maize Genomic Relationship Matrix (GRM)

This data set was included in Souza et al. (2017) and Cuevas et al
(2019) is from the Helix Seeds Company (HEL). Genotyping Array of 616 K
SNPs (Single Nucleotide Polymorphism) (Unterseer et al. 2014). Standard
quality controls (QC) were applied to the data, removing markers with a
Call Rate \$ 0.95. The remaining missing data in the lines were imputed
with the Synbreed package (Wimmer et al. 2015) using the algorithms from
the Beagle 4.0 software (Browning and Browning 2009). Markers with Minor
Allele Frequency (MAF) of \# 0.05 were removed. After applying QC,
52,811 SNPs were available to make the predictions. The phenotypic and
genomic data of inbred lines are credited to Helix Seeds Ltda. Company.

## Usage

``` r
data(maizeG)
```

## Format

An object of class `matrix` (inherits from `array`) with 150 rows and
150 columns.

## Examples

``` r
data(maizeG)
```
