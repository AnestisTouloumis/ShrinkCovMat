context("shrinkage towards the identity matrix")


set.seed(3)
p <- 20
n <- 10
datamat <- toeplitz(0.85 ^ seq(0, p - 1)) %*% matrix(rnorm(p * n), p, n)


test_that("shrinkage sphericity uncentered data", {
  sample_covariance_matrix <- cov(t(datamat))
  trace_sigma_hat <- sum(diag(sample_covariance_matrix))
  data_centered <- datamat - rowMeans(datamat)
  q <- sum(colSums(data_centered ^ 2) ^ 2) / (n - 1)
  trace_sigma_squared_hat <- (n - 1) / (n * (n - 2) * (n - 3)) *
    ((n - 1) * (n - 2) * sum(sample_covariance_matrix ^ 2) +
       trace_sigma_hat ^ 2 - n * q)
  lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
    (n * trace_sigma_squared_hat + (p - n + 1) / p * trace_sigma_hat ^ 2)
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(mean(diag(sample_covariance_matrix)), p)
  ans <- shrinkcovmat.equal(datamat)
  expect_equal(ans$Sigmahat, (1 - lambda_hat) * sample_covariance_matrix +
                 lambda_hat * target)
  expect_equal(ans$lambdahat, lambda_hat)
  expect_equal(ans$Sigmasample, sample_covariance_matrix)
  expect_equal(ans$Target, target)
})


test_that("shrinkage sphericity centered data", {
  sample_covariance_matrix <- tcrossprod(datamat) / n
  trace_sigma_hat <- sum(diag(sample_covariance_matrix))
  nu_hat <- trace_sigma_hat / p
  trace_sigma_squared_hat <- 0
  for (i in 1:(n - 1)) {
    trace_sigma_squared_hat <- sum(crossprod(datamat[, i],
                                             datamat[, (i + 1):n]) ^ 2) +
      trace_sigma_squared_hat
  }
  trace_sigma_squared_hat <- 2 * trace_sigma_squared_hat / n / (n - 1)
  lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat) /
    ((n + 1) * trace_sigma_squared_hat + (p - n) / p * trace_sigma_hat ^ 2)
  lambda_hat <- max(0, min(lambda_hat, 1))
  ans <- shrinkcovmat.equal(datamat, TRUE)
  target <- diag(nu_hat, p)
  expect_equal(ans$Sigmahat, (1 - lambda_hat) * sample_covariance_matrix +
                 lambda_hat * target)
  expect_equal(ans$lambdahat, lambda_hat)
  expect_equal(ans$Sigmasample, sample_covariance_matrix)
  expect_equal(ans$Target, target)
})

test_that("shrinkage sphericity centered argument", {
  expect_equal(shrinkcovmat.equal(datamat, "TRUE"),
               shrinkcovmat.equal(datamat, TRUE))
  expect_equal(shrinkcovmat.equal(datamat, "FALSE"),
               shrinkcovmat.equal(datamat, FALSE))
  expect_error(shrinkcovmat.equal(datamat, "iraklis"))
})

test_that("shrinkage sphericity sample size", {
  expect_error(shrinkcovmat.equal(datamat[, 1:3], FALSE))
  expect_error(shrinkcovmat.equal(datamat[, 1:2], FALSE))
  expect_error(shrinkcovmat.equal(datamat[, 1], FALSE))
  expect_error(shrinkcovmat.equal(datamat[, 1], TRUE))
})
