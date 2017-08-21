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
        DataCentered <- data - rowMeans(data)
        SigmaSample <- tcrossprod(DataCentered)/(N - 1)
        SigmaSampleVariances <- diag(SigmaSample)
        TraceSigmaHat <- sum(SigmaSampleVariances)
        Q <- sum(colSums(DataCentered^2)^2)/(N - 1)
        TraceSigmaSquaredHat <- (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 
            1) * (N - 2) * sum(SigmaSample^2) + (TraceSigmaHat)^2 - N * 
            Q)
        lambda1 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/
          (N * TraceSigmaSquaredHat + (p - N + 1)/p * TraceSigmaHat^2)
        lambda1 <- min(lambda1, 1)
        lambda2 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/
          (N * TraceSigmaSquaredHat + TraceSigmaHat^2 -
             2 * TraceSigmaHat * (N - 1) + p * (N - 1))
        lambda2 <- max(0, min(lambda2, 1))
        Sum1 <- Sum21 <- Sum22 <- Sum3 <- rep(0, p)
        for (i in 1:(N - 1)) {
            data2 <- matrix(data[, (i + 1):N], p, N - i)
            Sum1 <- rowSums(data[, i] * data2) + Sum1
            Sum21 <- rowSums(data[, i]^3 * data2) + Sum21
            Sum22 <- rowSums(data2^3 * data[, i]) + Sum22
            Sum3 <- rowSums(data[, i]^2 * data2^2) + Sum3
        }
        Term1 <- 2 * sum(Sum3)/N/(N - 1)
        Term2 <- 2 * (sum(Sum1 * rowSums(data^2)) - sum(Sum21 + Sum22))
        Term3 <- 4 * (sum(Sum1^2) - sum(Sum3) - Term2)
        Term2 <- Term2/N/(N - 1)/(N - 2)
        Term3 <- Term3/N/(N - 1)/(N - 2)/(N - 3)
        TraceDiagonalSigmaSquaredHat <- Term1 - 2 * Term2 + Term3
        lambda3 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat - 
                      2 * TraceDiagonalSigmaSquaredHat)/
          (N * TraceSigmaSquaredHat + TraceSigmaHat^2 -
             (N + 1) * TraceDiagonalSigmaSquaredHat)
        lambda3 <- max(0, min(lambda3, 1))
    } else {
        if (N < 2) 
            stop("The number of columns should be greater than 1")
        SigmaSample <- tcrossprod(data)/N
        SigmaSampleVariances <- diag(SigmaSample)
        TraceSigmaHat <- sum(SigmaSampleVariances)
        TraceSigmaSquaredHat <- 0
        TraceSigmaSquaredHat <- TraceDiagonalSigmaSquaredHat <- 0
        for (i in 1:(N - 1)) {
            TraceSigmaSquaredHat <- sum(crossprod(data[, i], data[, (i + 
                1):N])^2) + TraceSigmaSquaredHat
            TraceDiagonalSigmaSquaredHat <- sum((data[, i] * data[, (i + 
                1):N])^2) + TraceDiagonalSigmaSquaredHat
        }
        TraceSigmaSquaredHat <- 2 * TraceSigmaSquaredHat/N/(N - 1)
        TraceDiagonalSigmaSquaredHat <- 2 * TraceDiagonalSigmaSquaredHat/N/(N - 
            1)
        lambda1 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/((N + 1) * 
            TraceSigmaSquaredHat + (p - N)/p * TraceSigmaHat^2)
        lambda1 <- min(lambda1, 1)
        lambda2 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/((N + 1) * 
            TraceSigmaSquaredHat + TraceSigmaHat^2 - 2 * TraceSigmaHat * 
            N + p * N)
        lambda2 <- max(0, min(lambda2, 1))
        lambda3 <- (TraceSigmaHat^2 + TraceSigmaSquaredHat -
                      2 * TraceDiagonalSigmaSquaredHat)/
          ((N + 1) * TraceSigmaSquaredHat + TraceSigmaHat^2 -
             (N + 2) * TraceDiagonalSigmaSquaredHat)
        lambda3 <- max(0, min(lambda3, 1))
    }
    cat("OPTIMAL SHRINKAGE INTENSITIES FOR THE TARGET MATRIX WITH", "\n")
    cat("Equal variances   :", round(lambda1, 4), "\n")
    cat("Unit variances    :", round(lambda2, 4), "\n")
    cat("Unequal variances :", round(lambda3, 4), "\n")
    cat("\nSAMPLE VARIANCES", "\n")
    cat("Range   :", round(abs(diff(range(SigmaSampleVariances))), 4), "\n")
    cat("Average :", round(mean(SigmaSampleVariances), 4), "\n")
}
