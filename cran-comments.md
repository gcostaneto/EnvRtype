## Resubmission

This is a resubmission. In this version I have:

* Fixed the `kernel_model` help page, which showed an \arguments section
  without a \usage section. The documentation block was detached from the
  function definition, so roxygen2 did not emit \usage. The block is now
  attached and the regenerated Rd file contains both \usage and \arguments.
* Corrected example/documentation references that used the development
  package name instead of the final name (`EnvRtype`).

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* The NOTE reports possibly misspelled words in DESCRIPTION. These are
  domain terms in quantitative genetics and agronomy (envirotyping,
  enviromics, envirotypes, enviromes, covariable, evapotranspiration,
  bioclimatic, downscaled, agro) and the acronym FAO (Food and Agriculture
  Organization). All are spelled correctly.

## Test environments

* win-builder, R-release 4.6.1 (2026-06-24 ucrt): 1 NOTE
* win-builder, R-devel r90598 (2026-09-29 ucrt): 1 NOTE
* local Windows 11, R 4.6.0

## Notes

Some examples and vignette chunks use \donttest{} or eval = FALSE because
they require network access to NASA POWER, SoilGrids or WorldClim, or run
MCMC that would exceed reasonable check time.
