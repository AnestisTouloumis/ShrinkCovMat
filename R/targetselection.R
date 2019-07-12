#' Target Matrix Selection
#'
#' Implements the rule of thumb proposed by Touloumis (2015) for target matrix
#' selection. If the estimated optimal shrinkage intensities of the three
#' target matrices are of similar magnitude, then the average and the range of
#' the sample variances should be inspected in order to adopt the most
#' plausible target matrix.
#'
#' The rows of the data matrix \code{data} correspond to variables and the
#' columns to subjects.
#'
#' @param data a numeric matrix containing the data.
#' @param centered a logical indicating if the mean vector is the zero vector.
#' @return Prints the estimated optimal shrinkage intensities and the range and
#' the average of the sample variances.
#' @author Anestis Touloumis
#' @references Touloumis, A. (2015) Nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#' @examples
#' data(colon)
#' normal_group <- colon[, 1:40]
#' targetselection(normal_group)
#' ## Similar intensities, the range of the sample variances is small and the
#' ## average is not close to one. The scaled identity matrix seems to be the
#' ## most suitable target matrix for the normal group.
#'
#' tumor_group <- colon[, 41:62]
#' targetselection(tumor_group)
#' ## Similar intensities, the range of the sample variances is small and the
#' ## average is not close to one. The scaled identity matrix seems to be the
#' ## most suitable target matrix for the colon group.
#' @export
targetselection <- function(data, centered = FALSE) {
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  p <- nrow(data)
  n <- ncol(data)
  if (!centered) {
    if (n < 4) stop("The number of columns should be greater than 3")
    sample_variances <- apply(data, 1, var)
    trace_statistics <- trace_stats_uncentered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    trace_sigma_squared_hat <- trace_statistics[2]
    lambda_hat_shericity <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
      (n * trace_sigma_squared_hat + (p - n + 1) / p * trace_sigma_hat ^ 2)
    lambda_hat_shericity <- min(lambda_hat_shericity, 1)
    lambda_hat_identity <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
      (n * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
        2 * trace_sigma_hat * (n - 1) + p * (n - 1))
    lambda_hat_identity <- max(0, min(lambda_hat_identity, 1))
    trace_diagonal_sigma_sq_hat <- trace_statistics[3]
    lambda_hat_diagonal <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat -
      2 * trace_diagonal_sigma_sq_hat) /
      (n * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
        (n + 1) * trace_diagonal_sigma_sq_hat)
    lambda_hat_diagonal <- max(0, min(lambda_hat_diagonal, 1))
  } else {
    if (n < 2) stop("The number of columns should be greater than 1")
    sample_variances <- apply(data, 1, var)
    trace_statistics <- trace_stats_uncentered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    trace_sigma_squared_hat <- trace_statistics[2]
    trace_diagonal_sigma_sq_hat <- trace_statistics[3]
    lambda_hat_shericity <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
      ((n + 1) * trace_sigma_squared_hat + (p - n) / p * trace_sigma_hat ^ 2)
    lambda_hat_shericity <- min(lambda_hat_shericity, 1)
    lambda_hat_identity <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
      ((n + 1) * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
        2 * trace_sigma_hat * n + p * n)
    lambda_hat_identity <- max(0, min(lambda_hat_identity, 1))
    lambda_hat_diagonal <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat -
      2 * trace_diagonal_sigma_sq_hat) /
      ((n + 1) * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
        (n + 2) * trace_diagonal_sigma_sq_hat)
    lambda_hat_diagonal <- max(0, min(lambda_hat_diagonal, 1))
  }
  cat("OPTIMAL SHRInKAGE InTEnSITIES FOR THE TARGET MATRIX WITH", "\n")
  cat("Equal variances   :", round(lambda_hat_shericity, 4), "\n")
  cat("Unit variances    :", round(lambda_hat_identity, 4), "\n")
  cat("Unequal variances :", round(lambda_hat_diagonal, 4), "\n")
  cat("\nSAMPLE VARIAnCES", "\n")
  cat("Range   :", round(abs(diff(range(sample_variances))), 4), "\n")
  cat("Average :", round(mean(sample_variances), 4), "\n")
}
