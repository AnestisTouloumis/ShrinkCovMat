p <- sample(5:50, 1)
n <- sample(4:p, 1)
datamat <- matrix(rnorm(p * n), p, n)


test_that("sample covariance matrix - centered", {
  sample_covariance <- matrix(0, p, p)
  for (i in 1:n) {
    sample_covariance <- sample_covariance + tcrossprod(datamat[, i])
  }
  sample_covariance <- sample_covariance / n
  sample_covmat <- calculate_sample_covariance_matrix(datamat, TRUE, n)
  expect_equal(sample_covariance, sample_covmat)
})


test_that("sample covariance matrix - uncentered", {
  sample_covariance <- matrix(0, p, p)
  mean_vector <- rowMeans(datamat)
  for (i in 1:n) {
    sample_covariance <-
      sample_covariance + tcrossprod(datamat[, i] - mean_vector)
  }
  sample_covariance <- sample_covariance / (n - 1)
  sample_covmat <- calculate_sample_covariance_matrix(datamat, FALSE, n)
  expect_equal(sample_covariance, sample_covmat)
})


test_that("sample variances - centered", {
  sample_variances <- rep(0, p)
  for (i in 1:p) {
    for (j in 1:n) {
      sample_variances[i] <- sample_variances[i] + (datamat[i, j])^2
    }
  }
  sample_variances <- sample_variances / n
  sample_var <- calculate_sample_variances(datamat, TRUE)
  expect_equal(sample_variances, sample_var)
})


test_that("sample variances - uncentered", {
  sample_variances <- rep(0, p)
  mean_vector <- rowMeans(datamat)
  for (i in 1:p) {
    for (j in 1:n) {
      sample_variances[i] <-
        sample_variances[i] + (datamat[i, j] - mean_vector[i])^2
    }
  }
  sample_variances <- sample_variances / (n - 1)
  sample_var <- calculate_sample_variances(datamat, FALSE)
  expect_equal(sample_variances, sample_var)
})


test_that("calculate_lambda_hat - centered", {
  trace_stats <- calculate_trace_statistics(datamat, TRUE)
  sample_size <- n + 1
  lambda_diagonal <-
    calculate_lambda_hat(trace_stats, sample_size, p, "diagonal")
  lambda_spherical <-
    calculate_lambda_hat(trace_stats, sample_size, p, "spherical")
  lambda_identity <-
    calculate_lambda_hat(trace_stats, sample_size, p, "identity")
  lambdas <- targetselection(datamat, TRUE)
  expect_equal(lambda_identity, lambdas$lambda_hat_identity)
  expect_equal(lambda_diagonal, lambdas$lambda_hat_diagonal)
  expect_equal(lambda_spherical, lambdas$lambda_hat_spherical)
})


test_that("calculate_lambda_hat - uncentered", {
  trace_stats <- calculate_trace_statistics(datamat, FALSE)
  sample_size <- n
  lambda_diagonal <-
    calculate_lambda_hat(trace_stats, sample_size, p, "diagonal")
  lambda_spherical <-
    calculate_lambda_hat(trace_stats, sample_size, p, "spherical")
  lambda_identity <-
    calculate_lambda_hat(trace_stats, sample_size, p, "identity")
  lambdas <- targetselection(datamat, FALSE)
  expect_equal(lambda_identity, lambdas$lambda_hat_identity)
  expect_equal(lambda_diagonal, lambdas$lambda_hat_diagonal)
  expect_equal(lambda_spherical, lambdas$lambda_hat_spherical)
})


test_that("calculate_sigma_hat - centered", {
  sample_covariance <- calculate_sample_covariance_matrix(datamat, TRUE, n)
  sample_size <- n + 1
  trace_stats <- calculate_trace_statistics(datamat, TRUE)
  target_matrix <- calculate_target_matrix(datamat, TRUE, p, "diagonal")
  lambda_hat <- calculate_lambda_hat(trace_stats, sample_size, p, "diagonal")
  sigma_hat <- (1 - lambda_hat) * sample_covariance + lambda_hat * target_matrix
  sigma_matrix <-
    calculate_shrinkage_covariance_matrix(
      sample_covariance, lambda_hat, diag(target_matrix)
    )
  expect_equal(sigma_hat, sigma_matrix)
  target_matrix <- calculate_target_matrix(datamat, TRUE, p, "spherical")
  lambda_hat <- calculate_lambda_hat(trace_stats, sample_size, p, "spherical")
  sigma_hat <- (1 - lambda_hat) * sample_covariance + lambda_hat * target_matrix
  sigma_matrix <-
    calculate_shrinkage_covariance_matrix(
      sample_covariance, lambda_hat, diag(target_matrix)
    )
  expect_equal(sigma_hat, sigma_matrix)
  target_matrix <- calculate_target_matrix(datamat, TRUE, p, "identity")
  lambda_hat <- calculate_lambda_hat(trace_stats, sample_size, p, "identity")
  sigma_hat <- (1 - lambda_hat) * sample_covariance + lambda_hat * target_matrix
  sigma_matrix <-
    calculate_shrinkage_covariance_matrix(
      sample_covariance, lambda_hat, diag(target_matrix)
    )
  expect_equal(sigma_hat, sigma_matrix)
})


test_that("calculate_sigma_hat - uncentered", {
  sample_covariance <- calculate_sample_covariance_matrix(datamat, FALSE, n)
  sample_size <- n + 1
  trace_stats <- calculate_trace_statistics(datamat, FALSE)
  target_matrix <- calculate_target_matrix(datamat, FALSE, p, "diagonal")
  lambda_hat <- calculate_lambda_hat(trace_stats, sample_size, p, "diagonal")
  sigma_hat <- (1 - lambda_hat) * sample_covariance + lambda_hat * target_matrix
  sigma_matrix <-
    calculate_shrinkage_covariance_matrix(
      sample_covariance, lambda_hat, diag(target_matrix)
    )
  expect_equal(sigma_hat, sigma_matrix)
  target_matrix <- calculate_target_matrix(datamat, FALSE, p, "spherical")
  lambda_hat <- calculate_lambda_hat(trace_stats, sample_size, p, "spherical")
  sigma_hat <- (1 - lambda_hat) * sample_covariance + lambda_hat * target_matrix
  sigma_matrix <-
    calculate_shrinkage_covariance_matrix(
      sample_covariance, lambda_hat, diag(target_matrix)
    )
  expect_equal(sigma_hat, sigma_matrix)
  target_matrix <- calculate_target_matrix(datamat, FALSE, p, "identity")
  lambda_hat <- calculate_lambda_hat(trace_stats, sample_size, p, "identity")
  sigma_hat <- (1 - lambda_hat) * sample_covariance + lambda_hat * target_matrix
  sigma_matrix <-
    calculate_shrinkage_covariance_matrix(
      sample_covariance, lambda_hat, diag(target_matrix)
    )
  expect_equal(sigma_hat, sigma_matrix)
})
