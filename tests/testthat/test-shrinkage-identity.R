context("shrinkage towards the identity matrix")

set.seed(2)
p <- 20
N <- 5
datamat <- toeplitz(0.85^seq(0, p-1)) %*% matrix(rnorm(p*N), p, N)

test_that("uncentered data", {
  SampleCov <- cov(t(datamat))
  Y1N <- sum(diag(SampleCov))
  DataCentered <- datamat - rowMeans(datamat)
  Q <- sum(colSums(DataCentered^2)^2)/(N - 1)
  Y2N <- (N - 1)/(N * (N - 2) * (N - 3)) * 
    ((N - 1) * (N - 2) * sum(SampleCov^2) + (Y1N)^2 - N * Q)
  lambdahat <-  (Y1N^2 + Y2N)/ 
    (N * Y2N + Y1N^2 - 2 * Y1N * (N - 1) + p * (N - 1))
  lambdahat <- max(0, min(lambdahat, 1))
  Target <- diag(p)
  x <- shrinkcovmat.identity(datamat)
  expect_equal(x$Sigmahat, (1 - lambdahat) * SampleCov + lambdahat * Target)
  expect_equal(x$lambdahat, lambdahat)
  expect_equal(x$Sigmasample, SampleCov)
  expect_equal(x$Target, Target)
})


test_that("centered data", {
  SampleCov <- tcrossprod(datamat)/N
  Y1N <- sum(diag(SampleCov))
  Y2N <- 0
  for (i in 1:(N - 1)) {
    Y2N <- sum(crossprod(datamat[, i], datamat[, (i + 1):N])^2) + Y2N
  }
  Y2N <- 2 * Y2N/N/(N - 1)
  lambdahat <- (Y1N^2 + Y2N)/((N + 1) * Y2N + Y1N^2 - 2 * Y1N * N + p * N)
  lambdahat <- max(0, min(lambdahat, 1))
  Target <- diag(p)
  y <- shrinkcovmat.identity(datamat, centered = TRUE)
  expect_equal(y$Sigmahat, (1-lambdahat) * SampleCov + lambdahat * Target)
  expect_equal(y$lambdahat, lambdahat)
  expect_equal(y$Sigmasample, SampleCov)
  expect_equal(y$Target, Target)
})
