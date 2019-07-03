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
#' @references Touloumis, A. (2015) Nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#' @examples
#' data(colon)
#' NormalGroup <- colon[, 1:40]
#' TumorGroup <- colon[, 41:62]
#' Sigmahat.NormalGroup <- shrinkcovmat.equal(NormalGroup)
#' Sigmahat.NormalGroup
#' Sigmahat.TumorGroup <- shrinkcovmat.equal(TumorGroup)
#' Sigmahat.TumorGroup
#' @export
shrinkcovmat.equal <- function(data, centered = FALSE) {
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  p <- nrow(data)
  N <- ncol(data)
  centered <- as.logical(centered)
  if (centered != TRUE && centered != FALSE) {
    stop("'centered' must be either 'TRUE' or 'FALSE'")
  }
  if (!centered) {
    if (N < 4) stop("The number of columns should be greater than 3")
    sigma_sample <- cov(t(data))
    lambda_stats <- trace_stats_uncentered(data) # nolintr
    trace_sigma_hat <- lambda_stats[1]
    nu_hat <- trace_sigma_hat / p
    trace_sigma_squared_hat <- lambda_stats[2]
    lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
      (N * trace_sigma_squared_hat + (p - N + 1) /
        p * trace_sigma_hat ^ 2)
    lambda_hat <- min(lambda_hat, 1)
  } else {
    if (N < 2) stop("The number of columns should be greater than 1")
    sigma_sample <- tcrossprod(data) / N
    lambda_stats <- trace_stats_centered(data) # nolintr
    trace_sigma_hat <- lambda_stats[1]
    nu_hat <- trace_sigma_hat / p
    trace_sigma_squared_hat <- lambda_stats[2]
    lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
      ( (N + 1) * trace_sigma_squared_hat + (p - N) /
        p * trace_sigma_hat ^ 2)
    lambda_hat <- min(lambda_hat, 1)
  }
  if (lambda_hat < 1) {
    sigmahat <- (1 - lambda_hat) * sigma_sample +
      diag(nu_hat * lambda_hat, p)
  } else {
    sigmahat <- diag(lambda_hat * nu_hat, p)
  }
  target <- diag(nu_hat, p)
  ans <- list(
    Sigmahat = sigmahat, lambdahat = lambda_hat,
    Sigmasample = sigma_sample, Target = target,
    centered = centered
  )
  class(ans) <- "shrinkcovmathat"
  ans
}
