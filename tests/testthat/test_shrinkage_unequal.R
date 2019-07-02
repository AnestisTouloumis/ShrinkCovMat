context("shrinkage towards the diagonal matrix")

set.seed(3)
p <- 10
N <- 4
datamat <- toeplitz(0.85 ^ seq(0, p - 1)) %*% matrix(rnorm(p * N), p, N)


test_that("uncentered data", {
  sample_cov <- cov(t(datamat))
  data_centered <- datamat - rowMeans(datamat)
  sigma_sample_variances <- diag(sample_cov)
  Y1N <- sum(sigma_sample_variances)
  data_centered <- datamat - rowMeans(datamat)
  Q <- sum(colSums(data_centered ^ 2) ^ 2) / (N - 1)
  Y2N <- (N - 1) / (N * (N - 2) * (N - 3)) *
    ( (N - 1) * (N - 2) * sum(sample_cov ^ 2) + Y1N ^ 2 - N * Q)
  Sum1 <- Sum21 <- Sum22 <- Sum3 <- rep(0, p)
  for (i in 1:(N - 1)) {
    data2 <- matrix(datamat[, (i + 1):N], p, N - i)
    Sum1 <- rowSums(datamat[, i] * data2) + Sum1
    Sum21 <- rowSums(datamat[, i] ^ 3 * data2) + Sum21
    Sum22 <- rowSums(data2 ^ 3 * datamat[, i]) + Sum22
    Sum3 <- rowSums(datamat[, i] ^ 2 * data2 ^ 2) + Sum3
  }
  Term1 <- 2 * sum(Sum3) / N / (N - 1)
  Term2 <- 2 * (sum(Sum1 * rowSums(datamat ^ 2)) - sum(Sum21 + Sum22))
  Term3 <- 4 * (sum(Sum1 ^ 2) - sum(Sum3) - Term2)
  Term2 <- Term2 / N / (N - 1) / (N - 2)
  Term3 <- Term3 / N / (N - 1) / (N - 2) / (N - 3)
  Y3N <- Term1 - 2 * Term2 + Term3
  lambdahat <- (Y1N ^ 2 + Y2N - 2 * (1 - 1 / N) * Y3N) /
    (N * Y2N + Y1N ^ 2 - (N + 1 - 2 / N) * Y3N)
  lambdahat <- max(0, min(lambdahat, 1))
  x <- shrinkcovmat.unequal(datamat)
  target <- diag(sigma_sample_variances, p)
  expect_equal(x$Sigmahat, (1 - lambdahat) * sample_cov + lambdahat * target)
  expect_equal(x$lambdahat, lambdahat)
  expect_equal(x$Sigmasample, sample_cov)
  expect_equal(x$Target, target)
})



test_that("centered data", {
  sample_cov <- tcrossprod(datamat) / N
  sigma_sample_variances <- diag(sample_cov)
  Y1N <- sum(sigma_sample_variances)
  Y2N <- Y3N <- 0
  for (i in 1:(N - 1)) {
    Y2N <- sum(crossprod(datamat[, i], datamat[, (i + 1):N]) ^ 2) + Y2N
    Y3N <- sum( (datamat[, i] * datamat[, (i + 1):N]) ^ 2) + Y3N
  }
  Y2N <- 2 * Y2N / N / (N - 1)
  Y3N <- 2 * Y3N / N / (N - 1)
  lambdahat <- (Y1N ^ 2 + Y2N -  2 * (1 - 1 / (N + 1)) * Y3N) /
    ( (N + 1) * Y2N + Y1N ^ 2 - (N + 2 - 2 / (N + 1)) * Y3N)
  lambdahat <- max(0, min(lambdahat, 1))
  y <- shrinkcovmat.unequal(datamat, centered = TRUE)
  target <- diag(sigma_sample_variances, p)
  expect_equal(y$Sigmahat, (1 - lambdahat) * sample_cov + lambdahat * target)
  expect_equal(y$lambdahat, lambdahat)
  expect_equal(y$Sigmasample, sample_cov)
  expect_equal(y$Target, target)
})
