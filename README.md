<!-- README.md is generated from README.Rmd. Please edit that file -->
ShrinkCovMat: Estimation of a high-dimensional covariance matrix in R
=====================================================================

[![Travis-CI Build Status](https://travis-ci.org/AnestisTouloumis/ShrinkCovMat.svg?branch=master)](https://travis-ci.org/AnestisTouloumis/ShrinkCovMat) [![develVersion](https://img.shields.io/badge/devel%20version-1.1.3-brightgreen.svg?style=flat)](https://github.com/AnestisTouloumis/ShrinkCovMat) [![CRAN Version](http://www.r-pkg.org/badges/version/ShrinkCovMat?color=blue)](https://cran.r-project.org/package=ShrinkCovMat) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ShrinkCovMat?color=blue)](http://cranlogs.r-pkg.org/badges/grand-total/ShrinkCovMat)

Installation
------------

You can install the release version of the **ShrinkCovMat** package from CRAN:

``` r
install.packages("ShrinkCovMat")
```

of the development version from github:

``` r
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("AnestisTouloumis/ShrinkCovMat")
```

The source code for the **ShrinkCovMat** package is available on GitHub at

-   <https://github.com/AnestisTouloumis/ShrinkCovMat>.

Functions
---------

The **ShrinkCovMat** package provides three core function to estimate a covariance matrix.

-   `shrinkcovmat.equal` [![develVersion](https://img.shields.io/badge/devel%20version-1.1.3-blue.svg?style=flat)](https://github.com/AnestisTouloumis/ShrinkCovMat)
-   `shrinkcovmat.identity`
-   `shrinkcovmat.equal`

The utility function

-   `targetselection`

helps users to choose the target matrix.

To use these functions, first, you should load the package as follows.

``` r
library(ShrinkCovMat)
```
