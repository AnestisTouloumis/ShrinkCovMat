#' Shrinking the Sample Covariance Matrix Towards the Identity Matrix
#' 
#' Provides a nonparametric Stein-type shrinkage estimator of the covariance
#' matrix that is a linear combination of the sample covariance matrix and of
#' the identity matrix.
#' 
#' The rows of the data matrix \code{data} correspond to variables and the
#' columns to subjects.
#' 
#' @param data a numeric matrix containing the data.
#' @param centered a logical indicating if the mean vector is the zero vector.
#' @return Returns an object of the class 'shrinkcovmathat' that has
#' components: \item{Sigmahat}{The Stein-type shrinkage estimator of the
#' covariance matrix.} \item{lambdahat}{The estimated optimal shrinkage
#' intensity.} \item{Sigmasample}{The sample covariance matrix.}
#' \item{Target}{The target covariance matrix.} \item{centered}{If the data are
#' centered around their mean vector.}
#' @author Anestis Touloumis
#' @seealso \code{\link{shrinkcovmat.equal}} and
#' \code{\link{shrinkcovmat.unequal}}.
#' @references Touloumis, A. (2015) Nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#' @examples
#' data(colon)
#' NormalGroup <- colon[, 1:40]
#' TumorGroup <- colon[, 41:62]
#' Sigmahat.NormalGroup <- shrinkcovmat.identity(NormalGroup)
#' Sigmahat.NormalGroup
#' Sigmahat.TumorGroup <- shrinkcovmat.identity(TumorGroup)
#' Sigmahat.TumorGroup
#' @export
shrinkcovmat.identity <- function(data, centered = FALSE) {
    if (!is.matrix(data)) 
        data <- as.matrix(data)
    p <- nrow(data)
    N <- ncol(data)
    centered <- as.logical(centered)
    if (centered != TRUE && centered != FALSE) 
        stop("'centered' must be either 'TRUE' or 'FALSE'")
    if (!centered) {
        if (N < 4) 
            stop("The number of columns should be greater than 3")
        DataCentered <- data - rowMeans(data)
        SigmaSample <- tcrossprod(DataCentered)/(N - 1)
        TraceSigmaHat <- sum(diag(SigmaSample))
        Q <- sum(colSums(DataCentered^2)^2)/(N - 1)
        TraceSigmaSquaredHat <- (N - 1)/(N * (N - 2) * (N - 3)) * 
            ((N - 1) * (N - 2) * sum(SigmaSample^2) + (TraceSigmaHat)^2 - 
                N * Q)
        LambdaHat <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/(N * 
            TraceSigmaSquaredHat + TraceSigmaHat^2 - 2 * TraceSigmaHat * 
            (N - 1) + p * (N - 1))
        LambdaHat <- max(0, min(LambdaHat, 1))
    } else {
        if (N < 2) 
            stop("The number of columns should be greater than 1")
        SigmaSample <- tcrossprod(data)/N
        TraceSigmaHat <- sum(diag(SigmaSample))
        TraceSigmaSquaredHat <- 0
        for (i in 1:(N - 1)) TraceSigmaSquaredHat <- sum(crossprod(data[, 
            i], data[, (i + 1):N])^2) + TraceSigmaSquaredHat
        TraceSigmaSquaredHat <- 2 * TraceSigmaSquaredHat/N/(N - 1)
        LambdaHat <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/((N + 
            1) * TraceSigmaSquaredHat + TraceSigmaHat^2 - 2 * TraceSigmaHat * 
            N + p * N)
        LambdaHat <- max(0, min(LambdaHat, 1))
    }
    if (LambdaHat < 1) {
        SigmaHat <- (1 - LambdaHat) * SigmaSample + diag(LambdaHat, 
            p)
    } else SigmaHat <- diag(LambdaHat, p)
    Target <- diag(p)
    ans <- list(Sigmahat = SigmaHat, lambdahat = LambdaHat, Sigmasample = SigmaSample, 
        Target = Target, centered = centered)
    class(ans) <- "shrinkcovmathat"
    ans
}
70
70
70
