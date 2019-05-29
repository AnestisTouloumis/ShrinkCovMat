context("shrinkage towards the diagonal matrix")

set.seed(3)
p <- 10
N <- 4
datamat <- toeplitz(0.85^seq(0, p-1)) %*% matrix(rnorm(p*N), p, N)


test_that("uncentered data", {
  DataCentered <- datamat - rowMeans(datamat)
  SigmaSample <- tcrossprod(DataCentered)/(N - 1)
  SigmaSampleVariances <- diag(SigmaSample)
  TraceSigmaHat <- sum(SigmaSampleVariances)
  Q <- sum(colSums(DataCentered^2)^2)/(N - 1)
  TraceSigmaSquaredHat <- (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 1) * (N - 2) * sum(SigmaSample^2) + (TraceSigmaHat)^2 - N * Q)
  Sum1 <- Sum21 <- Sum22 <- Sum3 <- rep(0, p)
  for (i in 1:(N - 1)) {
    data2 <- matrix(datamat[, (i + 1):N], p, N - i)
    Sum1 <- rowSums(datamat[, i] * data2) + Sum1
    Sum21 <- rowSums(datamat[, i]^3 * data2) + Sum21
    Sum22 <- rowSums(data2^3 * datamat[, i]) + Sum22
    Sum3 <- rowSums(datamat[, i]^2 * data2^2) + Sum3
  }
  Term1 <- 2 * sum(Sum3)/N/(N - 1)
  Term2 <- 2 * (sum(Sum1 * rowSums(datamat^2)) - sum(Sum21 + Sum22))
  Term3 <- 4 * (sum(Sum1^2) - sum(Sum3) - Term2)
  Term2 <- Term2/N/(N - 1)/(N - 2)
  Term3 <- Term3/N/(N - 1)/(N - 2)/(N - 3)
  TraceDiagonalSigmaSquaredHat <- Term1 - 2 * Term2 + Term3
  lambda <- (TraceSigmaHat^2 + TraceSigmaSquaredHat - 2 * TraceDiagonalSigmaSquaredHat)/(N * TraceSigmaSquaredHat + TraceSigmaHat^2 - (N + 1) * TraceDiagonalSigmaSquaredHat)
  lambda <- max(0, min(lambda, 1))
  x <- shrinkcovmat.unequal(datamat)
  SampleCov <- cov(t(datamat))
  Target <- diag(diag(SampleCov), p)
  expect_equal(x$Sigmahat, (1-lambda) * SampleCov + lambda * Target)
  expect_equal(x$lambdahat, lambda)
  expect_equal(x$Sigmasample, SampleCov)
  expect_equal(x$Target, Target)
})



test_that("centered data", {
  SigmaSample <- tcrossprod(datamat)/N
  TraceSigmaHat <- sum(diag(SigmaSample))
  TraceSigmaSquaredHat <- TraceDiagonalSigmaSquaredHat <- 0
  for (i in 1:(N - 1)) {
    TraceSigmaSquaredHat <- sum(crossprod(datamat[, i], datamat[, (i + 
                                                               1):N])^2) + TraceSigmaSquaredHat
    TraceDiagonalSigmaSquaredHat <- sum((datamat[, i] * datamat[, (i + 
                                                               1):N])^2) + TraceDiagonalSigmaSquaredHat
  }
  TraceSigmaSquaredHat <- 2 * TraceSigmaSquaredHat/N/(N - 1)
  TraceDiagonalSigmaSquaredHat <- 2 * TraceDiagonalSigmaSquaredHat/(N*(N - 1))
  lambda <- (TraceSigmaHat^2 + TraceSigmaSquaredHat -  2 * TraceDiagonalSigmaSquaredHat)/ ((N + 1) * TraceSigmaSquaredHat + TraceSigmaHat^2 - (N + 2) * TraceDiagonalSigmaSquaredHat)
  lambda <- max(0, min(lambda, 1))
  y <- shrinkcovmat.unequal(datamat, centered = TRUE)
  Target <- diag(diag(SigmaSample), p)
  expect_equal(y$Sigmahat, (1-lambda) * SigmaSample + lambda * Target)
  expect_equal(y$lambdahat, lambda)
  expect_equal(y$Sigmasample, SigmaSample)
  expect_equal(y$Target, Target)
})