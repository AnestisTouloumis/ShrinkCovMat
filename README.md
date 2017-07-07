<!-- README.md is generated from README.Rmd. Please edit that file -->
ShrinkCovMat: Shrinkage Covariance Matrix Estimators
====================================================

[![Travis-CI Build Status](https://travis-ci.org/AnestisTouloumis/ShrinkCovMat.svg?branch=master)](https://travis-ci.org/AnestisTouloumis/ShrinkCovMat) [![Project Status: Active The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![CRAN Version](http://www.r-pkg.org/badges/version/ShrinkCovMat?color=blue)](https://CRAN.R-project.org/package=ShrinkCovMat) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ShrinkCovMat?color=blue)](http://cranlogs.r-pkg.org/badges/grand-total/ShrinkCovMat) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/ShrinkCovMat)](https://CRAN.R-project.org/package=ShrinkCovMat)

Author
------

Anestis Touloumis: <https://brighton.academia.edu/AnestisTouloumis>

School of Computing, Engineering and Mathematics, University of Brighton.

Installation
------------

You can install the release version of the **ShrinkCovMat**:

``` r
install.packages("ShrinkCovMat")
```

The source code for the release version of **ShrinkCovMat** is available on CRAN at:

-   <https://CRAN.R-project.org/package=ShrinkCovMat>

Or you can install the development version of **ShrinkCovMat**:

``` r
install.packages("devtools")  # if you have not installed 'devtools' package
devtools::install_github("AnestisTouloumis/ShrinkCovMat")
```

The source code for the development version of **ShrinkCovMat** is available on github at:

-   <https://github.com/AnestisTouloumis/ShrinkCovMat>

To use **ShrinkCovMat**, you should first load the package as follows:

``` r
library(ShrinkCovMat)
```

Usage
-----

This package provides nonparametric Stein-type shrinkage estimates of the covariance matrix that are suitable and statistically efficient when the number of variables is larger than the sample size. These estimators are non-singular and well-conditioned regardless of the dimensionality.

Each of the implemented shrinkage covariance matrix estimators is a convex linear combination of the sample covariance matrix and of a target matrix. Three options are considered for the target matrix:

-   the identity matrix (`shrinkcovmat.identity`),
-   the scaled identity matrix (`shrinkcovmat.equal`),
-   the diagonal matrix with diagonal elements the corresponding sample variances (`shrinkcovmat.unequal`).

Estimation of the corresponding optimal shrinkage intensities is discussed in Touloumis (2015). The utility function

-   `targetselection`

is designed to ease the selection of the target matrix.

Example
-------

Consider the colon cancer data example analyzed in Touloumis (2015). There are two groups: the normal group and the tumor group.

``` r
# load colon data
data(colon)
normal <- colon[, 1:40]
tumor <- colon[, 41:62]
```

We can use the `targetselection` function to chose between the three target matrices in each group. First, for the normal group

``` r
targetselection(normal)
#> OPTIMAL SHRINKAGE INTENSITIES FOR THE TARGET MATRIX WITH 
#> Equal variances   : 0.1401 
#> Unit variances    : 0.1125 
#> Unequal variances : 0.14 
#> 
#> SAMPLE VARIANCES 
#> Range   : 0.4714 
#> Average : 0.0882
```

we can see that the estimated shrinkage intensity for the scaled identity matrix is slightly larger than the other two and that the sample variances are of similar magnitude. Therefore we can choose the scaled identity matrix as the target matrix.

``` r
EstimatedCovarianceNormal <- shrinkcovmat.equal(normal)
EstimatedCovarianceNormal
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

Next, we follow a similar procedure for the tumor group

``` r
targetselection(tumor)
#> OPTIMAL SHRINKAGE INTENSITIES FOR THE TARGET MATRIX WITH 
#> Equal variances   : 0.1956 
#> Unit variances    : 0.1705 
#> Unequal variances : 0.1955 
#> 
#> SAMPLE VARIANCES 
#> Range   : 0.4226 
#> Average : 0.0958
```

Again, we should choose the scaled identity matrix as the target matrix since the estimated optimal shrinkage intensity for the scaled identity matrix is slightly larger than the other two and the sample variances are of similar magnitude.

``` r
EstimatedCovarianceTumor <- shrinkcovmat.equal(tumor)
EstimatedCovarianceTumor
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

How to cite
-----------

To cite *ShrinkCovMat* in publications, please use the following reference

    #> Touloumis, A. (2015) Nonparametric Stein-type Shrinkage Covariance Matrix Estimators in High-Dimensional Settings, Computational Statistics & Data Analysis 83, 251-261.

    #> @Article{,
    #>   title = {Nonparametric Stein-type Shrinkage Covariance Matrix Estimators in High-Dimensional Settings},
    #>   author = {{Anestis Touloumis}},
    #>   year = {2015},
    #>   journal = {Computational Statistics & Data Analysis},
    #>   volume = {83},
    #>   pages = {251--261},
    #> }

References
==========

Touloumis, Anestis. 2015. “Nonparametric Stein-Type Shrinkage Covariance Matrix Estimators in High-Dimensional Settings.” *Computational Statistics & Data Analysis* 83 (March): 251–61. doi:[10.1016/j.csda.2014.10.018](https://doi.org/10.1016/j.csda.2014.10.018).
