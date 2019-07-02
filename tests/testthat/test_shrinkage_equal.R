context("shrinkage towards the identity matrix")

set.seed(3)
p <- 20
N <- 10
datamat <- toeplitz(0.85 ^ seq(0, p - 1)) %*% matrix(rnorm(p * N), p, N)

test_that("uncentered data", {
  sample_cov <- cov(t(datamat))
  Y1N <- sum(diag(sample_cov))
  data_centered <- datamat - rowMeans(datamat)
  Q <- sum(colSums(data_centered ^ 2) ^ 2) / (N - 1)
  Y2N <- (N - 1) / (N * (N - 2) * (N - 3)) *
    ( (N - 1) * (N - 2) * sum(sample_cov ^ 2) + Y1N ^ 2 - N * Q)
  lambdahat <-  (Y1N ^ 2 + Y2N) / (N * Y2N + (p - N + 1) / p * Y1N ^ 2)
  lambdahat <- max(0, min(lambdahat, 1))
  target <- diag(mean(diag(sample_cov)), p)
  x <- shrinkcovmat.equal(datamat)
  expect_equal(x$Sigmahat, (1 - lambdahat) * sample_cov + lambdahat * target)
  expect_equal(x$lambdahat, lambdahat)
  expect_equal(x$Sigmasample, sample_cov)
  expect_equal(x$Target, target)
})


test_that("centered data", {
  sample_cov <- tcrossprod(datamat) / N
  Y1N <- sum(diag(sample_cov))
  nu_hat <- Y1N / p
  Y2N <- 0
  for (i in 1:(N - 1)) {
    Y2N <- sum(crossprod(datamat[, i], datamat[, (i + 1):N]) ^ 2) + Y2N
  }
  Y2N <- 2 * Y2N / N / (N - 1)
  lambdahat <- (Y1N ^ 2 + Y2N) / ( (N + 1) * Y2N + (p - N) / p * Y1N ^ 2)
  lambdahat <- max(0, min(lambdahat, 1))
  y <- shrinkcovmat.equal(datamat, TRUE)
  target <- diag(nu_hat, p)
  expect_equal(y$Sigmahat, (1 - lambdahat) * sample_cov + lambdahat * target)
  expect_equal(y$lambdahat, lambdahat)
  expect_equal(y$Sigmasample, sample_cov)
  expect_equal(y$Target, target)
})
