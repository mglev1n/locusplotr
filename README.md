
# locusplotr <img src='man/figures/logo.png' align="right" height="139" />
<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/mglev1n/locusplotr/workflows/R-CMD-check/badge.svg)](https://github.com/mglev1n/locusplotr/actions)
<!-- badges: end -->

The goal of `locusplotr` is to allow users to integrate genome-wide association study results with linkage disequlibirium data to create regional association plots surrounding a locus of interest. This package uses the `ld_extract_locuszoom` function to query the University of Michigan LocusZoom API (<https://portaldev.sph.umich.edu/>) to obtain linkage disequilibrium for a genetic variant and genomic region of interest. The `gg_locusplot` function allows the user to provide GWAS summary statistics for a region of interest and a reference variant, and will 1) use the `ld_extract_locuszoom` function to obtain LD information and 2) return a ggplot object with a regional association plot.

## Installation

You can install the released version of locusplotr from [GitHub](https://github.com/mglev1n/locusplotr) with:

``` r
devtools::install_github("mglev1n/locusplotr")
```
