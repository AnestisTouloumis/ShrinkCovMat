context("shrinkage towards the diagonal matrix")

set.seed(3)
p <- 10
n <- 4
datamat <- toeplitz(0.85 ^ seq(0, p - 1)) %*% matrix(rnorm(p * n), p, n)


test_that("uncentered data", {
  sample_cov <- cov(t(datamat))
  data_centered <- datamat - rowMeans(datamat)
  sigma_sample_variances <- diag(sample_cov)
  y_1n <- sum(sigma_sample_variances)
  data_centered <- datamat - rowMeans(datamat)
  q <- sum(colSums(data_centered ^ 2) ^ 2) / (n - 1)
  y_2n <- (n - 1) / (n * (n - 2) * (n - 3)) *
    ((n - 1) * (n - 2) * sum(sample_cov ^ 2) + y_1n ^ 2 - n * q)
  sum_1 <- sum_21 <- sum_22 <- sum_3 <- rep(0, p)
  for (i in 1:(n - 1)) {
    data2 <- matrix(datamat[, (i + 1):n], p, n - i)
    sum_1 <- rowSums(datamat[, i] * data2) + sum_1
    sum_21 <- rowSums(datamat[, i] ^ 3 * data2) + sum_21
    sum_22 <- rowSums(data2 ^ 3 * datamat[, i]) + sum_22
    sum_3 <- rowSums(datamat[, i] ^ 2 * data2 ^ 2) + sum_3
  }
  term_1 <- 2 * sum(sum_3) / n / (n - 1)
  term_2 <- 2 * (sum(sum_1 * rowSums(datamat ^ 2)) - sum(sum_21 + sum_22))
  term_3 <- 4 * (sum(sum_1 ^ 2) - sum(sum_3) - term_2)
  term_2 <- term_2 / n / (n - 1) / (n - 2)
  term_3 <- term_3 / n / (n - 1) / (n - 2) / (n - 3)
  y_3n <- term_1 - 2 * term_2 + term_3
  lambdahat <- (y_1n ^ 2 + y_2n - 2 * (1 - 1 / n) * y_3n) /
    (n * y_2n + y_1n ^ 2 - (n + 1 - 2 / n) * y_3n)
  lambdahat <- max(0, min(lambdahat, 1))
  x <- shrinkcovmat.unequal(datamat)
  target <- diag(sigma_sample_variances, p)
  expect_equal(x$Sigmahat, (1 - lambdahat) * sample_cov + lambdahat * target)
  expect_equal(x$lambdahat, lambdahat)
  expect_equal(x$Sigmasample, sample_cov)
  expect_equal(x$Target, target)
})



test_that("centered data", {
  sample_cov <- tcrossprod(datamat) / n
  sigma_sample_variances <- diag(sample_cov)
  y_1n <- sum(sigma_sample_variances)
  y_2n <- y_3n <- 0
  for (i in 1:(n - 1)) {
    y_2n <- sum(crossprod(datamat[, i], datamat[, (i + 1):n]) ^ 2) + y_2n
    y_3n <- sum((datamat[, i] * datamat[, (i + 1):n]) ^ 2) + y_3n
  }
  y_2n <- 2 * y_2n / n / (n - 1)
  y_3n <- 2 * y_3n / n / (n - 1)
  lambdahat <- (y_1n ^ 2 + y_2n - 2 * (1 - 1 / (n + 1)) * y_3n) /
    ((n + 1) * y_2n + y_1n ^ 2 - (n + 2 - 2 / (n + 1)) * y_3n)
  lambdahat <- max(0, min(lambdahat, 1))
  y <- shrinkcovmat.unequal(datamat, centered = TRUE)
  target <- diag(sigma_sample_variances, p)
  expect_equal(y$Sigmahat, (1 - lambdahat) * sample_cov + lambdahat * target)
  expect_equal(y$lambdahat, lambdahat)
  expect_equal(y$Sigmasample, sample_cov)
  expect_equal(y$Target, target)
})
