context("shrinkage towards the identity matrix")

set.seed(2)
p <- 20
n <- 5
datamat <- toeplitz(0.85 ^ seq(0, p - 1)) %*% matrix(rnorm(p * n), p, n)

test_that("shrinkage identity uncentered data", {
  sample_covariance_matrix <- cov(t(datamat))
  trace_sigma_hat <- sum(diag(sample_covariance_matrix))
  data_centered <- datamat - rowMeans(datamat)
  q <- sum(colSums(data_centered ^ 2) ^ 2) / (n - 1)
  trace_sigma_squared_hat <- (n - 1) / (n * (n - 2) * (n - 3)) *
    ((n - 1) * (n - 2) * sum(sample_covariance_matrix ^ 2) +
       trace_sigma_hat ^ 2 - n * q)
  lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
    (n * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
       2 * trace_sigma_hat * (n - 1) + p * (n - 1))
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(p)
  ans <- shrinkcovmat.identity(datamat)
  expect_equal(ans$Sigmahat, (1 - lambda_hat) * sample_covariance_matrix +
                 lambda_hat * target)
  expect_equal(ans$lambdahat, lambda_hat)
  expect_equal(ans$Sigmasample, sample_covariance_matrix)
  expect_equal(ans$Target, target)
})


test_that("shrinkage identity centered data", {
  sample_covariance_matrix <- tcrossprod(datamat) / n
  trace_sigma_hat <- sum(diag(sample_covariance_matrix))
  trace_sigma_squared_hat <- 0
  for (i in 1:(n - 1)) {
    trace_sigma_squared_hat <- sum(crossprod(datamat[, i],
                                             datamat[, (i + 1):n]) ^ 2) +
      trace_sigma_squared_hat
  }
  trace_sigma_squared_hat <- 2 * trace_sigma_squared_hat / n / (n - 1)
  lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
    ((n + 1) * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
       2 * trace_sigma_hat * n + p * n)
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(p)
  ans <- shrinkcovmat.identity(datamat, centered = TRUE)
  expect_equal(ans$Sigmahat, (1 - lambda_hat) * sample_covariance_matrix +
                 lambda_hat * target)
  expect_equal(ans$lambdahat, lambda_hat)
  expect_equal(ans$Sigmasample, sample_covariance_matrix)
  expect_equal(ans$Target, target)
})


test_that("shrinkage identity centered argument", {
  expect_equal(shrinkcovmat.identity(datamat, centered = "TRUE"),
               shrinkcovmat.identity(datamat, centered = TRUE))
  expect_equal(shrinkcovmat.identity(datamat, centered = "FALSE"),
               shrinkcovmat.identity(datamat, centered = FALSE))
  expect_error(shrinkcovmat.identity(datamat, centered = "iraklis"))
})


test_that("shrinkage identity sample size", {
  expect_error(shrinkcovmat.identity(datamat[, 1:3], FALSE))
  expect_error(shrinkcovmat.identity(datamat[, 1:2], FALSE))
  expect_error(shrinkcovmat.identity(datamat[, 1], FALSE))
  expect_error(shrinkcovmat.identity(datamat[, 1], TRUE))
})
