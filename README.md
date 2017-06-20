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

You can install the release version of the **ShrinkCovMat** package:

``` r
install.packages("ShrinkCovMat")
```

The source code for the release version of **ShrinkCovMat** is available on CRAN at:

-   <https://CRAN.R-project.org/package=ShrinkCovMat>

Or you can install the development version of the **ShrinkCovMat** package:

``` r
install.packages("devtools")  # if you have not installed 'devtools' package
devtools::install_github("AnestisTouloumis/ShrinkCovMat")
```

The source code for the development version of **ShrinkCovMat** is available on github at

-   <https://github.com/AnestisTouloumis/ShrinkCovMat>.

To use **ShrinkCovMat**, you should load the package as follows:

``` r
library(ShrinkCovMat)
```

Usage
-----

The **ShrinkCovMat** package provides nonparametric Stein-type shrinkage estimators of the covariance matrix that are suitable and statistically efficient when the number of variables is larger than the sample size. These estimators are non-singular and well-conditioned regardless of the dimensionality.

Each of the implemented shrinkage covariance matrix estimators is a convex linear combination of the sample covariance matrix and of a target matrix. Three options are considered for the target matrix:

-   the identity matrix (`shrinkcovmat.identity`),
-   the scaled identity matrix (`shrinkcovmat.equal`),
-   the diagonal matrix with diagonal elements the corresponding sample variances (`shrinkcovmat.unequal`).

Estimation of the corresponding optimal shrinkage intensities is discussed in (Touloumis 2015). The utility function

-   `targetselection`

is designed to ease the selection of the target matrix.

Example
-------

The following R code illustrates how to use the core functions of **ShrinkCovMat**.

``` r
data(colon)
normal.group <- colon[, 1:40]

Sigmahat1 <- shrinkcovmat.equal(normal.group)
Sigmahat1
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

Sigmahat2 <- shrinkcovmat.identity(normal.group)
Sigmahat2
#> SHRINKAGE ESTIMATION OF THE COVARIANCE MATRIX 
#> 
#> Estimated Optimal Shrinkage Intensity = 0.1125 
#> 
#> Estimated Covariance Matrix [1:5,1:5] =
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] 0.1406 0.0111 0.0105 0.0221 0.0181
#> [2,] 0.0111 0.1512 0.0379 0.0176 0.0041
#> [3,] 0.0105 0.0379 0.1512 0.0152 0.0046
#> [4,] 0.0221 0.0176 0.0152 0.1537 0.0094
#> [5,] 0.0181 0.0041 0.0046 0.0094 0.1495
#> 
#> Target Matrix [1:5,1:5] =
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1

Sigmahat3 <- shrinkcovmat.unequal(normal.group)
Sigmahat3
#> SHRINKAGE ESTIMATION OF THE COVARIANCE MATRIX 
#> 
#> Estimated Optimal Shrinkage Intensity = 0.14 
#> 
#> Estimated Covariance Matrix [1:5,1:5] =
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] 0.0044 0.0107 0.0101 0.0214 0.0175
#> [2,] 0.0107 0.0061 0.0368 0.0171 0.0040
#> [3,] 0.0101 0.0368 0.0061 0.0147 0.0045
#> [4,] 0.0214 0.0171 0.0147 0.0065 0.0091
#> [5,] 0.0175 0.0040 0.0045 0.0091 0.0058
#> 
#> Target Matrix [1:5,1:5] =
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] 0.0317 0.0000 0.0000 0.0000 0.0000
#> [2,] 0.0000 0.0437 0.0000 0.0000 0.0000
#> [3,] 0.0000 0.0000 0.0436 0.0000 0.0000
#> [4,] 0.0000 0.0000 0.0000 0.0465 0.0000
#> [5,] 0.0000 0.0000 0.0000 0.0000 0.0418
```

How to cite
-----------

``` r
citation("ShrinkCovMat")
#> 
#> To cite the R package 'ShrinkCovMat' in publications, please use:
#> 
#>   Touloumis, A. (2015) Nonparametric Stein-type Shrinkage
#>   Covariance Matrix Estimators in High-Dimensional Settings,
#>   Computational Statistics & Data Analysis 83, 251-261.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Nonparametric Stein-type Shrinkage Covariance Matrix Estimators in High-Dimensional Settings},
#>     author = {{Anestis Touloumis}},
#>     year = {2015},
#>     journal = {Computational Statistics & Data Analysis},
#>     volume = {83},
#>     pages = {251--261},
#>   }
```

References
==========

Touloumis, Anestis. 2015. “Nonparametric Stein-Type Shrinkage Covariance Matrix Estimators in High-Dimensional Settings.” *Computational Statistics & Data Analysis* 83 (March): 251–61. doi:[10.1016/j.csda.2014.10.018](https://doi.org/10.1016/j.csda.2014.10.018).
