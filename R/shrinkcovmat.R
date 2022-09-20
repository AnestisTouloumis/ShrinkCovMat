#' Linear Shrinkage of the Sample Covariance
#'
#' Provides a nonparametric Stein-type shrinkage estimator of the covariance
#' matrix that is a linear combination of the sample covariance matrix and of a
#' target matrix.
#'
#' Options for the target matrix include the \code{spherical} sample covariance
#' matrix (the diagonal matrix with diagonal elements the average of the sample
#' variances), the \code{diagonal} sample covariance matrix (the diagonal matrix
#' with diagonal elements the corresponding sample variances), and (c) the
#' \code{identity} matrix.
#'
#' The rows of the data matrix \code{data} correspond to variables and the
#' columns to subjects.
#'
#' @param data a numeric matrix containing the data.
#' @param centered a logical indicating if the mean vector is the zero vector.
#' @param target a character indicating the target matrix. Options include
#' 'spherical', 'identity' or 'diagonal'.
#' @return Returns an object of the class 'shrinkcovmathat' that has
#' components:
#' \item{Sigmahat}{The Stein-type shrinkage estimator of the
#' covariance matrix.}
#' \item{lambdahat}{The estimated optimal shrinkage
#' intensity.}
#' \item{Sigmasample}{The sample covariance matrix.}
#' \item{Target}{The target covariance matrix.}
#' \item{centered}{If the data are centered around their mean vector.}
#' @author Anestis Touloumis
#' @seealso \code{\link{targetselection}}.
#' @references Touloumis, A. (2015) nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#' @examples
#' data("colon")
#' normal_group <- colon[, 1:40]
#' tumor_group <- colon[, 41:62]
#' targetselection(normal_group)
#' sigma_hat_normal_group <- shrinkcovmat(normal_group)
#' sigma_hat_normal_group
#' targetselection(normal_group)
#' sigma_hat_tumor_group <- shrinkcovmat(tumor_group)
#' sigma_hat_tumor_group
#' @export
shrinkcovmat <- function(data, target = "spherical", centered = FALSE) { # nolint
  if (!is.matrix(data)) data <- as.matrix(data)
  p <- nrow(data)
  n <- ncol(data)
  targets <- c("spherical", "identity", "diagonal")
  if (!is.element(target, targets)) {
    stop("'target' must be 'spherical', 'identity' or 'diagonal'")
  }
  centered <- as.logical(centered)
  if (centered != TRUE && centered != FALSE) {
    stop("'centered' must be either 'TRUE' or 'FALSE'")
  }
  if (!centered) {
    if (n < 4) stop("The number of columns should be greater than 3")
  } else {
    if (n < 2) stop("The number of columns should be greater than 1")
  }
  if (target == "sphericity") {
    ans <- shrinkcovmat_equal(data, centered, p, n)
  } else if (target == "diagonal") {
    ans <- shrinkcovmat_unequal(data, centered, p, n)
  } else {
    ans <- shrinkcovmat_identity(data, centered, p, n)
  }
  ans
}
