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
        SigmaSample <- cov(t(data))
        lambda_stats <- optimal_intensities_uncentered(data, SigmaSample)
        TraceSigmaHat <-lambda_stats[1]
        TraceSigmaSquaredHat <- lambda_stats[2]
        TraceDiagonalSigmaSquaredHat <- lambda_stats[3]
        LambdaHat <- (TraceSigmaHat^2 + TraceSigmaSquaredHat - 2 * 
            TraceDiagonalSigmaSquaredHat)/(N * TraceSigmaSquaredHat + 
            TraceSigmaHat^2 - (N + 1) * TraceDiagonalSigmaSquaredHat)
        LambdaHat <- max(0, min(LambdaHat, 1))
    } else {
        if (N < 2) 
            stop("the number of columns should be greater than 1")
        SigmaSample <- tcrossprod(data)/N
        lambda_stats <- optimal_intensities_centered(data)
        TraceSigmaHat <-lambda_stats[1]
        TraceSigmaSquaredHat <- lambda_stats[2]
        TraceDiagonalSigmaSquaredHat <- lambda_stats[3]
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

