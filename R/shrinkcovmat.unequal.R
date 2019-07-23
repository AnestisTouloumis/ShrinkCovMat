#' Shrinking the Sample Covariance Matrix Towards a Diagonal Matrix with
#' Diagonal Elements the Sample Variances.
#'
#' Provides a nonparametric Stein-type shrinkage estimator of the covariance
#' matrix that is a linear combination of the sample covariance matrix and of
#' the diagonal matrix with elements the corresponding sample variances on the
#' diagonal and zeros elsewhere.
#'
#' The rows of the data matrix \code{data} correspond to variables and the
#' columns to subjects.
#'
#' @param data a numeric matrix containing the data.
#' @param centered a logical indicating if the vectors are centered around
#' their mean vector.
#' @return Returns an object of the class 'shrinkcovmathat' that has
#' components: \item{Sigmahat}{The Stein-type shrinkage estimator of the
#' covariance matrix.} \item{lambdahat}{The estimated optimal shrinkage
#' intensity.} \item{Sigmasample}{The sample covariance matrix.}
#' \item{Target}{The target covariance matrix.} \item{centered}{If the data are
#' centered around their mean vector.}
#' @author Anestis Touloumis
#' @seealso \code{\link{shrinkcovmat.equal}} and
#' \code{\link{shrinkcovmat.identity}}.
#' @references Touloumis, A. (2015) nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#' @examples
#' data(colon)
#' normal_group <- colon[, 1:40]
#' tumor_group <- colon[, 41:62]
#' sigma_hat_normal_group <- shrinkcovmat.unequal(normal_group)
#' sigma_hat_normal_group
#' sigma_hat_tumor_group <- shrinkcovmat.unequal(tumor_group)
#' sigma_hat_tumor_group
#' @export
shrinkcovmat.unequal <- function(data, centered = FALSE) { # nolint
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  p <- nrow(data)
  n <- ncol(data)
  centered <- as.logical(centered)
  if (centered != TRUE && centered != FALSE) {
    stop("'centered' must be either 'TRUE' or 'FALSE'")
  }
  if (!centered) {
    if (n < 4) stop("The number of columns should be greater than 3")
    sample_covariance_matrix <- cov(t(data))
    sample_variances <- apply(data, 1, var)
    trace_statistics <- trace_stats_uncentered(data) # nolint
    trace_sigma_hat <- trace_statistics[1]
    trace_sigma_squared_hat <- trace_statistics[2]
    trace_diagonal_sigma_sq_hat <- trace_statistics[3]
    lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat -
      (2 - 2 / n) * trace_diagonal_sigma_sq_hat) /
      (n * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
        (n + 1 - 2 / n) * trace_diagonal_sigma_sq_hat)
    lambda_hat <- max(0, min(lambda_hat, 1))
  } else {
    if (n < 2) stop("The number of columns should be greater than 1")
    sample_covariance_matrix <- tcrossprod(data) / n
    sample_variances <- apply(data, 1, function(x) mean(x ^ 2))
    trace_statistics <- trace_stats_centered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    trace_sigma_squared_hat <- trace_statistics[2]
    trace_diagonal_sigma_sq_hat <- trace_statistics[3]
    lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat -
      (2 - 2 / (n + 1)) * trace_diagonal_sigma_sq_hat) /
      ((n + 1) * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
        (n + 2 - 2 / (n + 1)) * trace_diagonal_sigma_sq_hat)
    lambda_hat <- max(0, min(lambda_hat, 1))
  }
  if (lambda_hat < 1) {
    sigma_hat <- (1 - lambda_hat) * sample_covariance_matrix +
      diag(lambda_hat * sample_variances, p)
  } else {
    sigma_hat <- diag(lambda_hat * sample_variances, p)
  }
  target <- diag(sample_variances, p)
  ans <- list(
    Sigmahat = sigma_hat, lambdahat = lambda_hat,
    Sigmasample = sample_covariance_matrix, Target = target,
    centered = centered
  )
  class(ans) <- "shrinkcovmathat"
  ans
}
