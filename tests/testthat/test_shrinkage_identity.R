p <- sample(5:50, 1)
n <- sample(4:p, 1)
datamat <- matrix(rnorm(p * n), p, n)

test_that("checking output with uncentered data", {
  sample_covariance_matrix <- cov(t(datamat))
  trace_sigma_hat <- sum(diag(sample_covariance_matrix))
  data_centered <- datamat - rowMeans(datamat)
  q <- sum(colSums(data_centered^2)^2) / (n - 1)
  trace_sigma_squared_hat <- (n - 1) / (n * (n - 2) * (n - 3)) *
    ((n - 1) * (n - 2) * sum(sample_covariance_matrix^2) +
      trace_sigma_hat^2 - n * q)
  lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
    (n * trace_sigma_squared_hat + trace_sigma_hat^2 -
      2 * trace_sigma_hat * (n - 1) + p * (n - 1))
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(p)
  ans <- shrinkcovmat(datamat, target = "identity", centered = FALSE)
  expect_equal(ans$Sigmahat, (1 - lambda_hat) * sample_covariance_matrix +
    lambda_hat * target)
  expect_equal(ans$lambdahat, lambda_hat)
  expect_equal(ans$Sigmasample, sample_covariance_matrix)
  expect_equal(ans$Target, target)
})


test_that("checking output with centered data", {
  sample_covariance_matrix <- tcrossprod(datamat) / n
  trace_sigma_hat <- sum(diag(sample_covariance_matrix))
  trace_sigma_squared_hat <- 0
  for (i in 1:(n - 1)) {
    trace_sigma_squared_hat <- sum(crossprod(
      datamat[, i],
      datamat[, (i + 1):n]
    )^2) +
      trace_sigma_squared_hat
  }
  trace_sigma_squared_hat <- 2 * trace_sigma_squared_hat / n / (n - 1)
  lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
    ((n + 1) * trace_sigma_squared_hat + trace_sigma_hat^2 -
      2 * trace_sigma_hat * n + p * n)
  lambda_hat <- max(0, min(lambda_hat, 1))
  target <- diag(p)
  ans <- shrinkcovmat(datamat, target = "identity", centered = TRUE)
  expect_equal(ans$Sigmahat, (1 - lambda_hat) * sample_covariance_matrix +
    lambda_hat * target)
  expect_equal(ans$lambdahat, lambda_hat)
  expect_equal(ans$Sigmasample, sample_covariance_matrix)
  expect_equal(ans$Target, target)
})


test_that("checking centered argument", {
  expect_equal(
    shrinkcovmat(datamat, target = "identity", centered = "TRUE"),
    shrinkcovmat(datamat, target = "identity", centered = TRUE)
  )
  expect_equal(
    shrinkcovmat(datamat, target = "identity", centered = "FALSE"),
    shrinkcovmat(datamat, target = "identity", centered = FALSE)
  )
  expect_error(shrinkcovmat(datamat, target = "identity", centered = "iraklis"))
})


test_that("checking sample size requirements", {
  expect_error(
    shrinkcovmat(datamat[, 1:3], target = "identity", centered = FALSE),
    "The number of columns should be greater than 3"
  )
  expect_error(
    shrinkcovmat(datamat[, 1:2], target = "identity", centered = FALSE),
    "The number of columns should be greater than 3"
  )
  expect_error(
    shrinkcovmat(datamat[, 1], target = "identity", centered = FALSE),
    "The number of columns should be greater than 3"
  )
  expect_error(
    shrinkcovmat(datamat[, 1], target = "identity", centered = TRUE),
    "The number of columns should be greater than 1"
  )
})
