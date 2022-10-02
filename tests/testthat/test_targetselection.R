p <- sample(5:50, 1)
n <- sample(4:p, 1)
datamat <- matrix(rnorm(p * n), p, n)


test_that("checking output with uncentered data", {
  optimal_hat_sphericity <- shrinkcovmat(datamat, target = "spherical", centered = FALSE)$lambdahat
  optimal_hat_identity <- shrinkcovmat(datamat, target = "identity", centered = FALSE)$lambdahat
  optimal_hat_diagonal <- shrinkcovmat(datamat, target = "diagonal", centered = FALSE)$lambdahat
  sample_variances <- apply(datamat, 1, var)
  range_variances <- diff(range(sample_variances))
  average_variances <- mean(sample_variances)
  select_target <- targetselection(datamat, FALSE)
  expect_equal(select_target$lambda_hat_spherical, optimal_hat_sphericity)
  expect_equal(select_target$lambda_hat_identity, optimal_hat_identity)
  expect_equal(select_target$lambda_hat_diagonal, optimal_hat_diagonal)
  expect_equal(select_target$range, range_variances)
  expect_equal(select_target$average, average_variances)
})


test_that("checking output with centered data", {
  optimal_hat_sphericity <- shrinkcovmat(datamat, target = "spherical", centered = TRUE)$lambdahat
  optimal_hat_identity <- shrinkcovmat(datamat, target = "identity", centered = TRUE)$lambdahat
  optimal_hat_diagonal <- shrinkcovmat(datamat, target = "diagonal", centered = TRUE)$lambdahat
  sample_variances <- apply(datamat, 1, function(x) mean(x^2))
  range_variances <- diff(range(sample_variances))
  average_variances <- mean(sample_variances)
  select_target <- targetselection(datamat, TRUE)
  expect_equal(select_target$lambda_hat_spherical, optimal_hat_sphericity)
  expect_equal(select_target$lambda_hat_identity, optimal_hat_identity)
  expect_equal(select_target$lambda_hat_diagonal, optimal_hat_diagonal)
  expect_equal(select_target$range, range_variances)
  expect_equal(select_target$average, average_variances)
})
