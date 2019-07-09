context("shrinkage towards the identity matrix")

set.seed(2)
p <- 20
n <- 5
datamat <- toeplitz(0.85 ^ seq(0, p - 1)) %*% matrix(rnorm(p * n), p, n)

test_that("uncentered data", {
  sample_cov <- cov(t(datamat))
  y_1n <- sum(diag(sample_cov))
  data_centered <- datamat - rowMeans(datamat)
  q <- sum(colSums(data_centered ^ 2) ^ 2) / (n - 1)
  y_2n <- (n - 1) / (n * (n - 2) * (n - 3)) *
    ((n - 1) * (n - 2) * sum(sample_cov ^ 2) + y_1n ^ 2 - n * q)
  lambda_hat <- (y_1n ^ 2 + y_2n) /
    (n * y_2n + y_1n ^ 2 - 2 * y_1n * (n - 1) + p * (n - 1))
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(p)
  x <- shrinkcovmat.identity(datamat)
  expect_equal(x$Sigmahat, (1 - lambda_hat) * sample_cov + lambda_hat * target)
  expect_equal(x$lambdahat, lambda_hat)
  expect_equal(x$Sigmasample, sample_cov)
  expect_equal(x$Target, target)
})


test_that("centered data", {
  sample_cov <- tcrossprod(datamat) / n
  y_1n <- sum(diag(sample_cov))
  y_2n <- 0
  for (i in 1:(n - 1)) {
    y_2n <- sum(crossprod(datamat[, i], datamat[, (i + 1):n]) ^ 2) + y_2n
  }
  y_2n <- 2 * y_2n / n / (n - 1)
  lambda_hat <- (y_1n ^ 2 + y_2n) /
    ((n + 1) * y_2n + y_1n ^ 2 - 2 * y_1n * n + p * n)
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(p)
  y <- shrinkcovmat.identity(datamat, centered = TRUE)
  expect_equal(y$Sigmahat, (1 - lambda_hat) * sample_cov + lambda_hat * target)
  expect_equal(y$lambdahat, lambda_hat)
  expect_equal(y$Sigmasample, sample_cov)
  expect_equal(y$Target, target)
})
