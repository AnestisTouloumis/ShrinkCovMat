p <- 10
n <- 5
datamat <- matrix(rnorm(p * n), p, n)

test_that("checking output with uncentered data", {
  sample_covariance_matrix <- cov(t(datamat))
  data_centered <- datamat - rowMeans(datamat)
  sigma_sample_variances <- diag(sample_covariance_matrix)
  trace_sigma_hat <- sum(sigma_sample_variances)
  data_centered <- datamat - rowMeans(datamat)
  q <- sum(colSums(data_centered ^ 2) ^ 2) / (n - 1)
  trace_sigma_squared_hat <- (n - 1) / (n * (n - 2) * (n - 3)) *
    ((n - 1) * (n - 2) * sum(sample_covariance_matrix ^ 2) +
       trace_sigma_hat ^ 2 - n * q)
  sum_1 <- sum_2 <- sum_3 <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) sum_1 <- sum_1 + sum(datamat[, i]^2 * datamat[, j]^2)
      for (k in 1:n) {
        if (i != j & i != k & j != k)
          sum_2 <- sum_2 + sum(datamat[, i]^2 * datamat[, j] * datamat[, k])
        for (l in 1:n) {
          if (i != j & i != k & i != l & j != k & j != l & k != l)
            sum_3 <- sum_3 + sum(datamat[, i] * datamat[, j] * datamat[, k] *
                                   datamat[, l])
        }
      }
    }
  }
  sum_1 <- sum_1 / n / (n - 1)
  sum_2 <- sum_2 / n / (n - 1) / (n - 2)
  sum_3 <- sum_3 / n / (n - 1) / (n - 2) / (n - 3)
  trace_diagonal_sigma_sq_hat <- sum_1 - 2 * sum_2 + sum_3
  lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat -
                   2 * (1 - 1 / n) * trace_diagonal_sigma_sq_hat) /
    (n * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
       (n + 1 - 2 / n) * trace_diagonal_sigma_sq_hat)
  lambda_hat <- max(0, min(lambda_hat, 1))
  ans <- shrinkcovmat.unequal(datamat)
  target <- diag(sigma_sample_variances, p)
  expect_equal(ans$Sigmahat, (1 - lambda_hat) * sample_covariance_matrix +
                 lambda_hat * target)
  expect_equal(ans$lambdahat, lambda_hat)
  expect_equal(ans$Sigmasample, sample_covariance_matrix)
  expect_equal(ans$Target, target)
})


test_that("checking output with centered data", {
  sample_covariance_matrix <- tcrossprod(datamat) / n
  sigma_sample_variances <- diag(sample_covariance_matrix)
  trace_sigma_hat <- sum(sigma_sample_variances)
  trace_sigma_squared_hat <- trace_diagonal_sigma_sq_hat <- 0
  for (i in 1:(n - 1)) {
    trace_sigma_squared_hat <- sum(crossprod(datamat[, i],
                                             datamat[, (i + 1):n]) ^ 2) +
      trace_sigma_squared_hat
    trace_diagonal_sigma_sq_hat <- sum((datamat[, i] *
                                          datamat[, (i + 1):n]) ^ 2) +
      trace_diagonal_sigma_sq_hat
  }
  trace_sigma_squared_hat <- 2 * trace_sigma_squared_hat / n / (n - 1)
  trace_diagonal_sigma_sq_hat <- 2 * trace_diagonal_sigma_sq_hat / n / (n - 1)
  lambda_hat <- (trace_sigma_hat ^ 2 + trace_sigma_squared_hat -
                   2 * (1 - 1 / (n + 1)) * trace_diagonal_sigma_sq_hat) /
    ((n + 1) * trace_sigma_squared_hat + trace_sigma_hat ^ 2 -
       (n + 2 - 2 / (n + 1)) * trace_diagonal_sigma_sq_hat)
  lambda_hat <- max(0, min(lambda_hat, 1))
  ans <- shrinkcovmat.unequal(datamat, centered = TRUE)
  target <- diag(sigma_sample_variances, p)
  expect_equal(ans$Sigmahat, (1 - lambda_hat) * sample_covariance_matrix +
                 lambda_hat * target)
  expect_equal(ans$lambdahat, lambda_hat)
  expect_equal(ans$Sigmasample, sample_covariance_matrix)
  expect_equal(ans$Target, target)
})


test_that("checking centered argument", {
  expect_equal(shrinkcovmat.unequal(datamat, "TRUE"),
               shrinkcovmat.unequal(datamat, TRUE))
  expect_equal(shrinkcovmat.unequal(datamat, "FALSE"),
               shrinkcovmat.unequal(datamat, FALSE))
  expect_error(shrinkcovmat.unequal(datamat, "iraklis"))
})


test_that("testing sample size requirements", {
  expect_error(shrinkcovmat.unequal(datamat[, 1:3], FALSE),
               "The number of columns should be greater than 3")
  expect_error(shrinkcovmat.unequal(datamat[, 1:2], FALSE),
               "The number of columns should be greater than 3")
  expect_error(shrinkcovmat.unequal(datamat[, 1], FALSE),
               "The number of columns should be greater than 3")
  expect_error(shrinkcovmat.unequal(datamat[, 1], TRUE),
               "The number of columns should be greater than 1")
})
