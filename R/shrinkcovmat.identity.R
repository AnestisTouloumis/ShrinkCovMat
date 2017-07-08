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
        TraceSigmaSquaredHat <- (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 
            1) * (N - 2) * sum(SigmaSample^2) + (TraceSigmaHat)^2 - N * 
            Q)
        LambdaHat <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/
          (N * TraceSigmaSquaredHat + 
             TraceSigmaHat^2 - 2 * TraceSigmaHat * (N - 1) + p * (N - 1))
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
        LambdaHat <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/((N + 1) * 
            TraceSigmaSquaredHat + TraceSigmaHat^2 - 2 * TraceSigmaHat * 
            N + p * N)
        LambdaHat <- max(0, min(LambdaHat, 1))
    }
    if (LambdaHat < 1) {
        SigmaHat <- (1 - LambdaHat) * SigmaSample
        diag(SigmaHat) <- LambdaHat + diag(SigmaHat)
    } else SigmaHat <- diag(LambdaHat, p)
    Target <- diag(p)
    ans <- list(Sigmahat = SigmaHat, lambdahat = LambdaHat, 
                Sigmasample = SigmaSample, Target = Target, centered = centered)
    class(ans) <- "shrinkcovmathat"
    ans
}
