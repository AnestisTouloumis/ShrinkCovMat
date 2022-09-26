p <- sample(5:50, 1)
n <- sample(4:p, 1)
datamat <- matrix(rnorm(p * n), p, n)

test_that("sample covariance matrix - centered", {
  sample_covariance <- matrix(0, p, p)
  for (i in 1:n) {
    sample_covariance <- sample_covariance + tcrossprod(datamat[, i])
    }
  sample_covariance <- sample_covariance/n
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
  sample_covariance <- sample_covariance/(n - 1)
  sample_covmat <- calculate_sample_covariance_matrix(datamat, FALSE, n)
  expect_equal(sample_covariance, sample_covmat)
})