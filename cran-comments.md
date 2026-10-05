## Resubmission

This is a resubmission of EnvRtype 0.1.0, which did not pass the incoming
pre-test. In this version I have:

* Exported `kernel_model()`. It was documented with \arguments but had no
  \usage section because it was not exported, which triggered the
  'Rd files without \usage' NOTE on Debian. It is a user-facing function
  and is now correctly exported.

## R CMD check results

Local check (Windows 11, R 4.6.0): 0 errors | 0 warnings | 0 notes

On CRAN's incoming checks a NOTE is expected for the new submission and
for possibly misspelled words in DESCRIPTION. Those words are domain terms
in quantitative genetics and agronomy (envirotyping, enviromics,
envirotypes, enviromes, covariable, evapotranspiration, bioclimatic,
downscaled, agro) and the acronym FAO (Food and Agriculture Organization).
All are spelled correctly.

## Test environments

* local Windows 11, R 4.6.0: 0 errors, 0 warnings, 0 notes
* win-builder, R-release 4.6.1 and R-devel: 1 NOTE (new submission,
  possibly misspelled words)
