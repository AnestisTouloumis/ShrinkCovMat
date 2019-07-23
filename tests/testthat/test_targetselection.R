p <- 20
n <- 10
datamat <- matrix(rnorm(p * n), p, n)


test_that("checking output with uncentered data", {
  optimal_hat_sphericity <- shrinkcovmat.equal(datamat)$lambdahat
  optimal_hat_identity <- shrinkcovmat.identity(datamat)$lambdahat
  optimal_hat_diagonal <- shrinkcovmat.unequal(datamat)$lambdahat
  sample_variances <- apply(datamat, 1, var)
  range_variances <- diff(range(sample_variances))
  average_variances <- mean(sample_variances)
  select_target <- targetselection(datamat, FALSE)
  expect_equal(select_target$optimal_sphericity, optimal_hat_sphericity)
  expect_equal(select_target$optimal_identity, optimal_hat_identity)
  expect_equal(select_target$optimal_diagonal, optimal_hat_diagonal)
  expect_equal(select_target$range, range_variances)
  expect_equal(select_target$average, average_variances)
})


test_that("checking output with centered data", {
  optimal_hat_sphericity <- shrinkcovmat.equal(datamat, TRUE)$lambdahat
  optimal_hat_identity <- shrinkcovmat.identity(datamat, TRUE)$lambdahat
  optimal_hat_diagonal <- shrinkcovmat.unequal(datamat, TRUE)$lambdahat
  sample_variances <- apply(datamat, 1, function(x) mean(x ^ 2))
  range_variances <- diff(range(sample_variances))
  average_variances <- mean(sample_variances)
  select_target <- targetselection(datamat, TRUE)
  expect_equal(select_target$optimal_sphericity, optimal_hat_sphericity)
  expect_equal(select_target$optimal_identity, optimal_hat_identity)
  expect_equal(select_target$optimal_diagonal, optimal_hat_diagonal)
  expect_equal(select_target$range, range_variances)
  expect_equal(select_target$average, average_variances)
})
