
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WWdPCR

<!-- badges: start -->
<!-- badges: end -->

WWdPCR is a package to analyse dPCR duplex assays. The main use cases
are: computing confidence intervals for the variant proprotion in a dPCR
duplex assay, fitting a logistic growth model to the variant proportion
with time series duplex assay dPCR data (also including background
prevalence) to estimate growth rates. The logistic growth functions
(that include background prevalence) can also be applied to binomial
data.

If you use WWdPCR, please cite Caduff et al. 2021
[(https://doi.org/10.1101/2021.08.22.21262024)](https://doi.org/10.1101/2021.08.22.21262024).

## Installation

WWdPCR is not yet on CRAN.

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cbg-ethz/WWdPCR")
```

## Examples

To see examples on how to use WWdPCR, compile the vignettes on
installing:

``` r
# install.packages("devtools")
# install.packages("rmarkdown")
devtools::install_github("cbg-ethz/WWdPCR", build_vignettes=TRUE)
```

And then run the vignettes:

``` r
vignette("fitting_dPCR_data")
vignette("fitting_binomial_data")
vignette("example_confints")
```
