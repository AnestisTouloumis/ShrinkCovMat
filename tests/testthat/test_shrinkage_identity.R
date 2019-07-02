context("shrinkage towards the identity matrix")

set.seed(2)
p <- 20
N <- 5
datamat <- toeplitz(0.85 ^ seq(0, p - 1)) %*% matrix(rnorm(p * N), p, N)

test_that("uncentered data", {
  sample_cov <- cov(t(datamat))
  Y1N <- sum(diag(sample_cov))
  data_centered <- datamat - rowMeans(datamat)
  Q <- sum(colSums(data_centered ^ 2) ^ 2) / (N - 1)
  Y2N <- (N - 1) / (N * (N - 2) * (N - 3)) *
    ( (N - 1) * (N - 2) * sum(sample_cov ^ 2) + Y1N ^ 2 - N * Q)
  lambda_hat <-  (Y1N ^ 2 + Y2N) /
    (N * Y2N + Y1N ^ 2 - 2 * Y1N * (N - 1) + p * (N - 1))
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(p)
  x <- shrinkcovmat.identity(datamat)
  expect_equal(x$Sigmahat, (1 - lambda_hat) * sample_cov + lambda_hat * target)
  expect_equal(x$lambdahat, lambda_hat)
  expect_equal(x$Sigmasample, sample_cov)
  expect_equal(x$Target, target)
})


test_that("centered data", {
  sample_cov <- tcrossprod(datamat) / N
  Y1N <- sum(diag(sample_cov))
  Y2N <- 0
  for (i in 1:(N - 1)) {
    Y2N <- sum(crossprod(datamat[, i], datamat[, (i + 1):N]) ^ 2) + Y2N
  }
  Y2N <- 2 * Y2N / N / (N - 1)
  lambda_hat <- (Y1N ^ 2 + Y2N) /
    ( (N + 1) * Y2N + Y1N ^ 2 - 2 * Y1N * N + p * N)
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(p)
  y <- shrinkcovmat.identity(datamat, centered = TRUE)
  expect_equal(y$Sigmahat, (1 - lambda_hat) * sample_cov + lambda_hat * target)
  expect_equal(y$lambdahat, lambda_hat)
  expect_equal(y$Sigmasample, sample_cov)
  expect_equal(y$Target, target)
})
