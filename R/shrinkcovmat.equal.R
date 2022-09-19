#' Shrinking the Sample Covariance Matrix Towards a Sphericity Matrix
#'
#' Provides a nonparametric Stein-type shrinkage estimator of the covariance
#' matrix that is a linear combination of the sample covariance matrix and of a
#' diagonal matrix with the average of the sample variances on the diagonal and
#' zeros elsewhere.
#'
#' The rows of the data matrix \code{data} correspond to variables and the
#' columns to subjects.
#'
#' @param data a numeric matrix containing the data.
#' @param centered a logical indicating if the mean vector is the zero vector.
#' @return Returns an object of the class 'shrinkcovmathat' that has
#' components: \item{SigmaHat}{The Stein-type shrinkage estimator of the
#' covariance matrix.} \item{lambdahat}{The estimated optimal shrinkage
#' intensity.} \item{sigmasample}{The sample covariance matrix.}
#' \item{Target}{The target covariance matrix.} \item{centered}{If the data are
#' centered around their mean vector.}
#' @author Anestis Touloumis
#' @seealso \code{\link{shrinkcovmat.unequal}} and
#' \code{\link{shrinkcovmat.identity}}.
#' @references Touloumis, A. (2015) nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#'
#' @name shrinkcovmat.equal-deprecated
#' @usage shrinkcovmat.equal(data, centered)
#' @seealso \code{\link{ShrinkCovMat-deprecated}}
#' @keywords internal
NULL
#' @rdname ShrinkCovMat-deprecated
#' @section \code{shrinkcovmat.equal}:
#' For \code{shrinkcovmat.equal}, use \code{\link{shrinkcovmat}}.
#'
#' @export
shrinkcovmat.equal <- function(data, centered = FALSE) { # nolint
  .Deprecated(msg = "Use instead function 'shrinkcovmat' with argument 'target' equal to 'sphericity'")
  if (!is.matrix(data)) data <- as.matrix(data)
  p <- nrow(data)
  n <- ncol(data)
  centered <- as.logical(centered)
  if (centered != TRUE && centered != FALSE) {
    stop("'centered' must be either 'TRUE' or 'FALSE'")
  }
  if (!centered) {
    if (n < 4) stop("The number of columns should be greater than 3")
    sample_covariance_matrix <- cov(t(data))
    trace_statistics <- trace_stats_uncentered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    nu_hat <- trace_sigma_hat / p
    trace_sigma_squared_hat <- trace_statistics[2]
    lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
      (n * trace_sigma_squared_hat + (p - n + 1) / p * trace_sigma_hat^2)
    lambda_hat <- min(lambda_hat, 1)
  } else {
    if (n < 2) stop("The number of columns should be greater than 1")
    sample_covariance_matrix <- tcrossprod(data) / n
    trace_statistics <- trace_stats_centered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    nu_hat <- trace_sigma_hat / p
    trace_sigma_squared_hat <- trace_statistics[2]
    lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
      ((n + 1) * trace_sigma_squared_hat + (p - n) / p * trace_sigma_hat^2)
    lambda_hat <- min(lambda_hat, 1)
  }
  if (lambda_hat < 1) {
    sigmahat <- (1 - lambda_hat) * sample_covariance_matrix +
      diag(nu_hat * lambda_hat, p)
  } else {
    sigmahat <- diag(lambda_hat * nu_hat, p)
  }
  target <- diag(nu_hat, p)
  ans <- list(
    Sigmahat = sigmahat, lambdahat = lambda_hat,
    Sigmasample = sample_covariance_matrix, Target = target,
    centered = centered
  )
  class(ans) <- "shrinkcovmathat"
  ans
}
