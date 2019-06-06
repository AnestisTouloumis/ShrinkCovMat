context("shrinkage towards the diagonal matrix")

set.seed(3)
p <- 10
N <- 4
datamat <- toeplitz(0.85^seq(0, p-1)) %*% matrix(rnorm(p*N), p, N)


test_that("uncentered data", {
  SampleCov <- cov(t(datamat))
  DataCentered <- datamat - rowMeans(datamat)
  SigmaSampleVariances <- diag(SampleCov)
  Y1N <- sum(SigmaSampleVariances)
  DataCentered <- datamat - rowMeans(datamat)
  Q <- sum(colSums(DataCentered^2)^2)/(N - 1)
  Y2N <- (N - 1)/(N * (N - 2) * (N - 3)) * 
    ((N - 1) * (N - 2) * sum(SampleCov^2) + (Y1N)^2 - N * Q)
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
  Y3N <- Term1 - 2 * Term2 + Term3
  lambdahat <- (Y1N^2 + Y2N - 2 * Y3N)/(N * Y2N + Y1N^2 - (N + 1) * Y3N)
  lambdahat <- max(0, min(lambdahat, 1))
  x <- shrinkcovmat.unequal(datamat)
  Target <- diag(SigmaSampleVariances, p)
  expect_equal(x$Sigmahat, (1-lambdahat) * SampleCov + lambdahat * Target)
  expect_equal(x$lambdahat, lambdahat)
  expect_equal(x$Sigmasample, SampleCov)
  expect_equal(x$Target, Target)
})



test_that("centered data", {
  SigmaCov <- tcrossprod(datamat)/N
  SigmaSampleVariances <- diag(SigmaCov)
  Y1N <- sum(SigmaSampleVariances)
  Y2N <- Y3N <- 0
  for (i in 1:(N - 1)) {
    Y2N <- sum(crossprod(datamat[, i], datamat[, (i + 1):N])^2) + Y2N
    Y3N <- sum((datamat[, i] * datamat[, (i + 1):N])^2) + Y3N
  }
  Y2N <- 2 * Y2N/N/(N - 1)
  Y3N <- 2 * Y3N/N/(N - 1)
  lambdahat <- (Y1N^2 + Y2N -  2 * Y3N)/ ((N + 1) * Y2N + Y1N^2 - (N + 2) * Y3N)
  lambdahat <- max(0, min(lambdahat, 1))
  y <- shrinkcovmat.unequal(datamat, centered = TRUE)
  Target <- diag(SigmaSampleVariances, p)
  expect_equal(y$Sigmahat, (1 - lambdahat) * SigmaCov + lambdahat * Target)
  expect_equal(y$lambdahat, lambdahat)
  expect_equal(y$Sigmasample, SigmaCov)
  expect_equal(y$Target, Target)
})
