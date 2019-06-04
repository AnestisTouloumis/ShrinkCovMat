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
#' NormalGroup <- colon[, 1:40]
#' targetselection(NormalGroup)
#' ## Similar intensities, the range of the sample variances is small and the average
#' ## is not close to one. The scaled identity matrix seems to be the most suitable
#' ## target matrix for the normal group
#' 
#' TumorGroup <- colon[, 41:62]
#' targetselection(TumorGroup)
#' ## Similar intensities, the range of the sample variances is small and the average
#' ## is not close to one. The scaled identity matrix seems to be the most suitable
#' ## target matrix for the colon group
#' @export
targetselection <- function(data, centered = FALSE) {
    if (!is.matrix(data)) 
        data <- as.matrix(data)
    p <- nrow(data)
    N <- ncol(data)
    if (!centered) {
        if (N < 4) 
            stop("The number of columns should be greater than 3")
        SigmaSample <- cov(t(data))
        SigmaSampleVariances <- apply(data, 1, var)
        lambda_stats <- optimal_intensities_uncentered(data, SigmaSample)
        TraceSigmaHat <- lambda_stats[1]
        TraceSigmaSquaredHat <- lambda_stats[2]
        lambda1 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/(N * 
            TraceSigmaSquaredHat + (p - N + 1)/p * TraceSigmaHat^2)
        lambda1 <- min(lambda1, 1)
        lambda2 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/(N * 
            TraceSigmaSquaredHat + TraceSigmaHat^2 - 2 * TraceSigmaHat * 
            (N - 1) + p * (N - 1))
        lambda2 <- max(0, min(lambda2, 1))
        TraceDiagonalSigmaSquaredHat <- lambda_stats[3]
        lambda3 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat - 2 * 
            TraceDiagonalSigmaSquaredHat)/(N * TraceSigmaSquaredHat + 
            TraceSigmaHat^2 - (N + 1) * TraceDiagonalSigmaSquaredHat)
        lambda3 <- max(0, min(lambda3, 1))
    } else {
        if (N < 2) 
            stop("The number of columns should be greater than 1")
        SigmaSampleVariances <- apply(data, 1, var)
        lambda_stats <- optimal_intensities_centered(data)
        TraceSigmaHat <- lambda_stats[1]
        TraceSigmaSquaredHat <- lambda_stats[2]
        TraceDiagonalSigmaSquaredHat <- lambda_stats[3]
        lambda1 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/((N + 
            1) * TraceSigmaSquaredHat + (p - N)/p * TraceSigmaHat^2)
        lambda1 <- min(lambda1, 1)
        lambda2 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/((N + 
            1) * TraceSigmaSquaredHat + TraceSigmaHat^2 - 2 * TraceSigmaHat * 
            N + p * N)
        lambda2 <- max(0, min(lambda2, 1))
        lambda3 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat - 2 * 
            TraceDiagonalSigmaSquaredHat)/((N + 1) * TraceSigmaSquaredHat + 
            TraceSigmaHat^2 - (N + 2) * TraceDiagonalSigmaSquaredHat)
        lambda3 <- max(0, min(lambda3, 1))
    }
    cat("OPTIMAL SHRINKAGE INTENSITIES FOR THE TARGET MATRIX WITH", 
        "\n")
    cat("Equal variances   :", round(lambda1, 4), "\n")
    cat("Unit variances    :", round(lambda2, 4), "\n")
    cat("Unequal variances :", round(lambda3, 4), "\n")
    cat("\nSAMPLE VARIANCES", "\n")
    cat("Range   :", round(abs(diff(range(SigmaSampleVariances))), 
        4), "\n")
    cat("Average :", round(mean(SigmaSampleVariances), 4), "\n")
}
