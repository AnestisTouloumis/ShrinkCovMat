#' Shrinkage Covariance Matrix Estimators
#'
#' Provides nonparametric Stein-type shrinkage estimators of the covariance
#' matrix that are suitable and statistically efficient when the number of
#' variables is larger than the sample size. These estimators are non-singular
#' and well-conditioned regardless of the dimensionality.
#'
#' Each of the implemented shrinkage covariance matrix estimators is a convex
#' linear combination of the sample covariance matrix and of a target matrix.
#'
#' The function \code{\link{shrinkcovmat}} implements three options for the
#' target matrix: (a) spherical sample covariance matrix, i.e. the diagonal
#' matrix with diagonal elements the average of the sample variances, (b)
#' diagonal sample covariance matrix, i.e. the diagonal matrix with diagonal
#' elements the corresponding sample variances, and (c) the identity matrix
#' (\code{identity}). The optimal shrinkage intensity determines how much the
#' sample covariance matrix will be shrunk towards the selected target matrix.
#'
#' Estimation of the corresponding optimal shrinkage intensities is discussed
#' in \cite{Touloumis (2015)}. The function \code{\link{targetselection}} is
#' designed to ease the selection of the target matrix.
#'
#' @name ShrinkCovMat-package
#' @aliases ShrinkCovMat-package ShrinkCovMat
#' @docType package
#' @author Anestis Touloumis
#'
#' Maintainer: Anestis Touloumis <A.Touloumis@@brighton.ac.uk>
#' @references Touloumis, A. (2015) Nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#' @useDynLib ShrinkCovMat, .registration = TRUE
#' @import Rcpp
#' @importFrom stats cov var
#' @keywords package
"_PACKAGE"
# > [1] '_PACKAGE'
