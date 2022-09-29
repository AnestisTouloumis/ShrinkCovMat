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
#' @return Prints the estimated optimal shrinkage intensities, the range and
#' average of the sample variances and returns an object of the class
#' 'targetsel' that has components: \item{optimal_sphericity}{The estimated
#' optimal intensity for a target matrix with equal variances.}
#' \item{optimal_identity}{The estimated optimal shrinkage
#' intensity for the identity target matrix.} \item{optimal_diagonal}{The
#' estimated optimal intensity for a target matrix with unequal variances.}
#' \item{range}{The range of the sample variances.} \item{average}{The average
#' of the sample variances.}
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
  if (!is.matrix(data)) data <- as.matrix(data)
  p <- nrow(data)
  n <- ncol(data)
  centered <- as.logical(centered)
  if (centered != TRUE && centered != FALSE) {
    stop("'centered' must be either 'TRUE' or 'FALSE'")
  }
  if (!centered) {
    if (n < 4) stop("The number of columns should be greater than 3")
  } else {
    if (n < 2) stop("The number of columns should be greater than 1")
  }
  sample_variances <- calculate_sample_variances(data, centered)
  trace_statistics <- calculate_trace_statistics(data, centered)
  sample_size <- ifelse(centered, n + 1, n)
  lambda_hat_sphericity <- 
    calculate_lambda_hat(
      trace_statistics, sample_size, p, "spherical"
      )
  lambda_hat_identity <-
    calculate_lambda_hat(
      trace_statistics, sample_size, p, "identity"
      )
  lambda_hat_diagonal <-
    calculate_lambda_hat(
      trace_statistics, sample_size, p, "diagonal"
    ) 
  ans <- list(
    optimal_sphericity = lambda_hat_sphericity,
    optimal_identity = lambda_hat_identity,
    optimal_diagonal = lambda_hat_diagonal,
    range = abs(diff(range(sample_variances))),
    average = mean(sample_variances)
  )
  class(ans) <- "targetsel"
  ans
}
