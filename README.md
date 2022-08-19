
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HERON

<!-- badges: start -->
<!-- badges: end -->

The goal of HERON (**H**ierarchical **E**pitope P**R**Otein Bi**N**ding)
is to anayzing peptide binding array data.

## Installation

You can install the released version of HERON from
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

<img src="man/figures/README-flowchart_probe_pvalues-1.png" width="100%" />

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(HERON)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
