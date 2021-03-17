
<!-- README.md is generated from README.Rmd. Please edit that file -->
# meeqtl

<!-- badges: start -->
<!-- badges: end -->
The goal of meeqtl is to ...

## Installation

You can install the released version of meeqtl from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("meeqtl")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yzeng-lol/meeqtl")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(meeqtl)
#> Loading required package: tidyverse
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
#> ✓ tibble  3.1.0     ✓ dplyr   1.0.5
#> ✓ tidyr   1.1.3     ✓ stringr 1.4.0
#> ✓ readr   1.4.0     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
#> Loading required package: foreach
#> 
#> Attaching package: 'foreach'
#> The following objects are masked from 'package:purrr':
#> 
#>     accumulate, when
#> Warning: replacing previous import 'IRanges::colnames' by 'MatrixEQTL::colnames'
#> when loading 'meeqtl'
#> Warning: replacing previous import 'IRanges::rownames' by 'MatrixEQTL::rownames'
#> when loading 'meeqtl'
#> Warning: replacing previous import 'IRanges::collapse' by 'dplyr::collapse' when
#> loading 'meeqtl'
#> Warning: replacing previous import 'IRanges::union' by 'dplyr::union' when
#> loading 'meeqtl'
#> Warning: replacing previous import 'IRanges::slice' by 'dplyr::slice' when
#> loading 'meeqtl'
#> Warning: replacing previous import 'IRanges::intersect' by 'dplyr::intersect'
#> when loading 'meeqtl'
#> Warning: replacing previous import 'IRanges::setdiff' by 'dplyr::setdiff' when
#> loading 'meeqtl'
#> Warning: replacing previous import 'IRanges::desc' by 'dplyr::desc' when loading
#> 'meeqtl'
#> 
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

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

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
