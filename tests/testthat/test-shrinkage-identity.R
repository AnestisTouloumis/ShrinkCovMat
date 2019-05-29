context("shrinkage towards the identity matrix")

set.seed(2)
p <- 20
N <- 5
datamat <- toeplitz(0.85^seq(0, p-1)) %*% matrix(rnorm(p*N), p, N)

test_that("uncentered data", {
  Target <- diag(p)
  SampleCov <- cov(t(datamat))
  DataCentered <- datamat - rowMeans(datamat)
  SigmaSample <- tcrossprod(DataCentered)/(N - 1)
  SigmaSampleVariances <- diag(SigmaSample)
  TraceSigmaHat <- sum(SigmaSampleVariances)
  Q <- sum(colSums(DataCentered^2)^2)/(N - 1)
  TraceSigmaSquaredHat <- (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 1) * (N - 2) * sum(SigmaSample^2) + (TraceSigmaHat)^2 - N * Q)
  lambda <-  (TraceSigmaHat^2 + TraceSigmaSquaredHat)/ (N * TraceSigmaSquaredHat + TraceSigmaHat^2 - 2 * TraceSigmaHat * (N - 1) + p * (N - 1))
  lambda <- max(0, min(lambda, 1))
  x <- shrinkcovmat.identity(datamat)
  expect_equal(x$Sigmahat, (1-lambda) * SampleCov + lambda * Target)
  expect_equal(x$lambdahat, lambda)
  expect_equal(x$Sigmasample, SampleCov)
  expect_equal(x$Target, Target)
})


test_that("centered data", {
  SigmaSample <- tcrossprod(datamat)/N
  TraceSigmaHat <- sum(diag(SigmaSample))
  TraceSigmaSquaredHat <- 0
  for (i in 1:(N - 1)) TraceSigmaSquaredHat <- sum(crossprod(datamat[, i], datamat[, (i + 1):N])^2) + TraceSigmaSquaredHat
  TraceSigmaSquaredHat <- 2 * TraceSigmaSquaredHat/N/(N - 1)
  lambda <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/((N + 1) * TraceSigmaSquaredHat + TraceSigmaHat^2 - 2 * TraceSigmaHat * N + p * N)
  lambda <- max(0, min(lambda, 1))
  Target <- diag(p)
  SampleCov <- tcrossprod(datamat)/N
  y <- shrinkcovmat.identity(datamat, centered = TRUE)
  expect_equal(y$Sigmahat, (1-lambda) * SampleCov + lambda * Target)
  expect_equal(y$lambdahat, lambda)
  expect_equal(y$Sigmasample, SampleCov)
  expect_equal(y$Target, Target)
})