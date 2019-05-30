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
#' @references Touloumis, A. (2015) Nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#' @examples
#' data(colon)
#' NormalGroup <- colon[, 1:40]
#' TumorGroup <- colon[, 41:62]
#' Sigmahat.NormalGroup <- shrinkcovmat.unequal(NormalGroup)
#' Sigmahat.NormalGroup
#' Sigmahat.TumorGroup <- shrinkcovmat.unequal(TumorGroup)
#' Sigmahat.TumorGroup
#' @export 
shrinkcovmat.unequal <- function(data, centered = FALSE) {
    if (!is.matrix(data)) 
        data <- as.matrix(data)
    p <- nrow(data)
    N <- ncol(data)
    centered <- as.logical(centered)
    if (centered != TRUE && centered != FALSE) 
        stop("'centered' must be either 'TRUE' or 'FALSE'")
    if (!centered) {
        if (N < 4) 
            stop("the number of columns should be greater than 3")
        DataCentered <- data - rowMeans(data)
        SigmaSample <- tcrossprod(DataCentered)/(N - 1)
        TraceSigmaHat <- sum(diag(SigmaSample))
        Q <- sum(colSums(DataCentered^2)^2)/(N - 1)
        TraceSigmaSquaredHat <- (N - 1)/(N * (N - 2) * (N - 3)) * 
            ((N - 1) * (N - 2) * sum(SigmaSample^2) + (TraceSigmaHat)^2 - 
                N * Q)
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
        LambdaHat <- (TraceSigmaHat^2 + TraceSigmaSquaredHat - 2 * 
            TraceDiagonalSigmaSquaredHat)/(N * TraceSigmaSquaredHat + 
            TraceSigmaHat^2 - (N + 1) * TraceDiagonalSigmaSquaredHat)
        LambdaHat <- max(0, min(LambdaHat, 1))
    } else {
        if (N < 2) 
            stop("the number of columns should be greater than 1")
        SigmaSample <- tcrossprod(data)/N
        TraceSigmaHat <- sum(diag(SigmaSample))
        TraceSigmaSquaredHat <- TraceDiagonalSigmaSquaredHat <- 0
        for (i in 1:(N - 1)) {
            TraceSigmaSquaredHat <- sum(crossprod(data[, i], data[, 
                (i + 1):N])^2) + TraceSigmaSquaredHat
            TraceDiagonalSigmaSquaredHat <- sum((data[, i] * data[, 
                (i + 1):N])^2) + TraceDiagonalSigmaSquaredHat
        }
        TraceSigmaSquaredHat <- 2 * TraceSigmaSquaredHat/N/(N - 1)
        TraceDiagonalSigmaSquaredHat <- 2 * TraceDiagonalSigmaSquaredHat/(N * 
            (N - 1))
        LambdaHat <- (TraceSigmaHat^2 + TraceSigmaSquaredHat - 2 * 
            TraceDiagonalSigmaSquaredHat)/((N + 1) * TraceSigmaSquaredHat + 
            TraceSigmaHat^2 - (N + 2) * TraceDiagonalSigmaSquaredHat)
        LambdaHat <- max(0, min(LambdaHat, 1))
    }
    DiagonalSigmaSample <- diag(SigmaSample)
    if (LambdaHat < 1) {
        SigmaHat <- (1 - LambdaHat) * SigmaSample + diag(LambdaHat * 
            DiagonalSigmaSample, p)
    } else SigmaHat <- diag(LambdaHat * DiagonalSigmaSample, p)
    Target <- diag(DiagonalSigmaSample, p)
    ans <- list(Sigmahat = SigmaHat, lambdahat = LambdaHat, Sigmasample = SigmaSample, 
        Target = Target, centered = centered)
    class(ans) <- "shrinkcovmathat"
    ans
}
70
70
70
