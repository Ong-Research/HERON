
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HERON

<!-- badges: start -->
<!-- badges: end -->

The goal of HERON (**H**ierarchical **E**pitope p**RO**tein bi**N**ding)
is to analyze peptide binding array data measured from nimblegen.

## Installation

You can install the released version of HERON from
[bioconductor](https://www.bioconductor.org/) with:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("HERON")
```

Or you can install the released version of HERON from
[github](https://github.com/Ong-Research/HERON) with:

``` r
devtools::install_github("Ong-Research/HERON")
```

And the development version from
[GitHub](https://github.com/Ong-Research/HERON) with:

``` r
# install.packages("devtools")
devtools::install_github("Ong-Research/HERON")
```
