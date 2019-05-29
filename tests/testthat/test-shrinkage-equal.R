context("shrinkage towards the identity matrix")

set.seed(3)
p <- 20
N <- 10
datamat <- toeplitz(0.85^seq(0, p-1)) %*% matrix(rnorm(p*N), p, N)

test_that("uncentered data", {
  DataCentered <- datamat - rowMeans(datamat)
  SigmaSample <- tcrossprod(DataCentered)/(N - 1)
  SigmaSampleVariances <- diag(SigmaSample)
  TraceSigmaHat <- sum(SigmaSampleVariances)
  Q <- sum(colSums(DataCentered^2)^2)/(N - 1)
  TraceSigmaSquaredHat <- (N - 1)/(N * (N - 2) * (N - 3)) * ((N - 1) * (N - 2) * sum(SigmaSample^2) + (TraceSigmaHat)^2 - N * Q)
  lambda <-  (TraceSigmaHat^2 + TraceSigmaSquaredHat)/(N * TraceSigmaSquaredHat + (p - N + 1)/p * TraceSigmaHat^2)
  lambda <- max(0, min(lambda, 1))
  SampleCov <- cov(t(datamat))
  Target <- diag(mean(diag(SampleCov)), p)
  x <- shrinkcovmat.equal(datamat)
  expect_equal(x$Sigmahat, (1-lambda) * SampleCov + lambda * Target)
  expect_equal(x$lambdahat, lambda)
  expect_equal(x$Sigmasample, SampleCov)
  expect_equal(x$Target, Target)
})


test_that("centered data", {
  SigmaSample <- tcrossprod(datamat)/N
  TraceSigmaHat <- sum(diag(SigmaSample))
  NuHat <- TraceSigmaHat/p
  TraceSigmaSquaredHat <- 0
  for (i in 1:(N - 1)) TraceSigmaSquaredHat <- sum(crossprod(datamat[, i], datamat[, (i + 1):N])^2) + TraceSigmaSquaredHat
  TraceSigmaSquaredHat <- 2 * TraceSigmaSquaredHat/N/(N - 1)
  lambda <- (TraceSigmaHat^2 + TraceSigmaSquaredHat)/((N + 1) * TraceSigmaSquaredHat + (p - N)/p * TraceSigmaHat^2)
  lambda <- min(lambda, 1)
  y <- shrinkcovmat.equal(datamat, TRUE)
  SampleCov <- tcrossprod(datamat)/N
  Target <- diag(mean(diag(SampleCov)), p)
  expect_equal(y$Sigmahat, (1-lambda) * SampleCov + lambda * Target)
  expect_equal(y$lambdahat, lambda)
  expect_equal(y$Sigmasample, SampleCov)
  expect_equal(y$Target, Target)
})
