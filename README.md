
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ShrinkCovMat: Shrinkage Covariance Matrix Estimators

[![Github
version](https://img.shields.io/badge/GitHub%20-1.4.3-green.svg)](%22commits/master%22)
[![Travis-CI Build
Status](https://travis-ci.org/AnestisTouloumis/ShrinkCovMat.svg?branch=master)](https://travis-ci.org/AnestisTouloumis/ShrinkCovMat)
[![Project Status: Active The project has reached a stable, usable state
and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/AnestisTouloumis/ShrinkCovMat/branch/master/graph/badge.svg)](https://codecov.io/gh/AnestisTouloumis/ShrinkCovMat?branch=master)

[![CRAN
Version](https://www.r-pkg.org/badges/version/ShrinkCovMat?color=blue)](https://CRAN.R-project.org/package=ShrinkCovMat)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/ShrinkCovMat?color=blue)](https://cranlogs.r-pkg.org/badges/grand-total/ShrinkCovMat)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/ShrinkCovMat)](https://CRAN.R-project.org/package=ShrinkCovMat)

## Installation

You can install the release version of `ShrinkCovMat`:

``` r
install.packages("ShrinkCovMat")
```

The source code for the release version of `ShrinkCovMat` is available
on CRAN at:

  - <https://CRAN.R-project.org/package=ShrinkCovMat>

Or you can install the development version of `ShrinkCovMat`:

``` r
# install.packages('devtools')
devtools::install_github("AnestisTouloumis/ShrinkCovMat")
```

The source code for the development version of `ShrinkCovMat` is
available on github at:

  - <https://github.com/AnestisTouloumis/ShrinkCovMat>

To use `ShrinkCovMat`, you should first load the package as follows:

``` r
library(ShrinkCovMat)
```

## Usage

This package provides the nonparametric Stein-type shrinkage covariance
matrix estimators proposed by Touloumis (2015). These estimators are
suitable and statistically efficient regardless of the dimensionality.

Each of the implemented shrinkage covariance matrix estimators is a
convex linear combination of the sample covariance matrix and of a
target matrix. Three options are considered for the target matrix:

  - the identity matrix (`shrinkcovmat.identity`),
  - the scaled identity matrix (`shrinkcovmat.equal`),
  - the diagonal matrix with diagonal elements the corresponding sample
    variances (`shrinkcovmat.unequal`).

Estimation of the corresponding optimal shrinkage intensities is
discussed in Touloumis (2015).

The utility function `targetselection` is designed to ease the selection
of the target matrix.

## Example

Consider the colon cancer data example analyzed in Touloumis (2015). The
data consists of two tissue groups: the normal tissue group and the
tumor tissue group.

``` r
data(colon)
normal_group <- colon[, 1:40]
tumor_group <- colon[, 41:62]
```

To decide the target matrix for covariance matrix of the normal group,
inspect the output of the function `targetselection`:

``` r
targetselection(normal_group)
#> OPTIMAL SHRINKAGE INTENSITIES FOR THE TARGET MATRIX WITH 
#> Equal variances   : 0.1401 
#> Unit variances    : 0.1125 
#> Unequal variances : 0.14 
#> 
#> SAMPLE VARIANCES 
#> Range   : 0.4714 
#> Average : 0.0882
```

The estimated optimal shrinkage intensity for the scaled identity matrix
is slightly larger than the other two and the sample variances appear to
be of similar magnitude. Thus, the scaled identity matrix seems to be
the most appropriate target for the covariance estimated. The resulting
covariance matrix estimator is:

``` r
estimated_covariance_normal <- shrinkcovmat.equal(normal_group)
estimated_covariance_normal
#> SHRINKAGE ESTIMATION OF THE COVARIANCE MATRIX 
#> 
#> Estimated Optimal Shrinkage Intensity = 0.1401 
#> 
#> Estimated Covariance Matrix [1:5,1:5] =
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] 0.0396 0.0107 0.0101 0.0214 0.0175
#> [2,] 0.0107 0.0499 0.0368 0.0171 0.0040
#> [3,] 0.0101 0.0368 0.0499 0.0147 0.0045
#> [4,] 0.0214 0.0171 0.0147 0.0523 0.0091
#> [5,] 0.0175 0.0040 0.0045 0.0091 0.0483
#> 
#> Target Matrix [1:5,1:5] =
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] 0.0882 0.0000 0.0000 0.0000 0.0000
#> [2,] 0.0000 0.0882 0.0000 0.0000 0.0000
#> [3,] 0.0000 0.0000 0.0882 0.0000 0.0000
#> [4,] 0.0000 0.0000 0.0000 0.0882 0.0000
#> [5,] 0.0000 0.0000 0.0000 0.0000 0.0882
```

We follow a similar procedure for the tumor group:

``` r
targetselection(tumor_group)
#> OPTIMAL SHRINKAGE INTENSITIES FOR THE TARGET MATRIX WITH 
#> Equal variances   : 0.1956 
#> Unit variances    : 0.1705 
#> Unequal variances : 0.1955 
#> 
#> SAMPLE VARIANCES 
#> Range   : 0.4226 
#> Average : 0.0958
```

As before, we may conclude that the scaled identity matrix seems to be
the most suitable target matrix. The resulting covariance matrix
estimate for the tumor group is:

``` r
estimated_covariance_tumor <- shrinkcovmat.equal(tumor_group)
estimated_covariance_tumor
#> SHRINKAGE ESTIMATION OF THE COVARIANCE MATRIX 
#> 
#> Estimated Optimal Shrinkage Intensity = 0.1956 
#> 
#> Estimated Covariance Matrix [1:5,1:5] =
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] 0.0490 0.0179 0.0170 0.0195 0.0052
#> [2,] 0.0179 0.0450 0.0265 0.0092 0.0034
#> [3,] 0.0170 0.0265 0.0465 0.0084 0.0031
#> [4,] 0.0195 0.0092 0.0084 0.0498 0.0036
#> [5,] 0.0052 0.0034 0.0031 0.0036 0.0361
#> 
#> Target Matrix [1:5,1:5] =
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] 0.0958 0.0000 0.0000 0.0000 0.0000
#> [2,] 0.0000 0.0958 0.0000 0.0000 0.0000
#> [3,] 0.0000 0.0000 0.0958 0.0000 0.0000
#> [4,] 0.0000 0.0000 0.0000 0.0958 0.0000
#> [5,] 0.0000 0.0000 0.0000 0.0000 0.0958
```

## How to cite

``` 

To cite the R package 'ShrinkCovMat' in publications, please use:

  Touloumis, A. (2015) Nonparametric Stein-type Shrinkage Covariance
  Matrix Estimators in High-Dimensional Settings, Computational
  Statistics & Data Analysis 83, 251-261.

A BibTeX entry for LaTeX users is

  @Article{,
    title = {Nonparametric Stein-type Shrinkage Covariance Matrix Estimators in High-Dimensional Settings},
    author = {{Anestis Touloumis}},
    year = {2015},
    journal = {Computational Statistics & Data Analysis},
    volume = {83},
    pages = {251--261},
  }
```

# References

<div id="refs" class="references">

<div id="ref-Touloumis2015">

Touloumis, A. (2015) Nonparametric Stein-type Shrinkage Covariance
Matrix Estimators in High-Dimensional Settings. *Computational
Statistics & Data Analysis*, **83**, 251–261.

</div>

</div>
