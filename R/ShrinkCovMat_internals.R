shrinkcovmat_unequal <- function(data, centered, p, n) {
  if (!centered) {
    if (n < 4) stop("The number of columns should be greater than 3")
    sample_covariance_matrix <- cov(t(data))
    sample_variances <- apply(data, 1, var)
    trace_statistics <- trace_stats_uncentered(data) # nolint
    trace_sigma_hat <- trace_statistics[1]
    trace_sigma_squared_hat <- trace_statistics[2]
    trace_diagonal_sigma_sq_hat <- trace_statistics[3]
    lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat -
                     (2 - 2 / n) * trace_diagonal_sigma_sq_hat) /
      (n * trace_sigma_squared_hat + trace_sigma_hat^2 -
         (n + 1 - 2 / n) * trace_diagonal_sigma_sq_hat)
    lambda_hat <- max(0, min(lambda_hat, 1)) 
    } else {
    if (n < 2) stop("The number of columns should be greater than 1")
    sample_covariance_matrix <- tcrossprod(data) / n
    sample_variances <- apply(data, 1, function(x) mean(x^2))
    trace_statistics <- trace_stats_centered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    trace_sigma_squared_hat <- trace_statistics[2]
    trace_diagonal_sigma_sq_hat <- trace_statistics[3]
    lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat -
                     (2 - 2 / (n + 1)) * trace_diagonal_sigma_sq_hat) /
      ((n + 1) * trace_sigma_squared_hat + trace_sigma_hat^2 -
         (n + 2 - 2 / (n + 1)) * trace_diagonal_sigma_sq_hat)
    lambda_hat <- max(0, min(lambda_hat, 1))
    }
  if (lambda_hat < 1) {
    sigma_hat <- (1 - lambda_hat) * sample_covariance_matrix +
      diag(lambda_hat * sample_variances, p)
  } else {
    sigma_hat <- diag(lambda_hat * sample_variances, p)
  }
  target <- diag(sample_variances, p)
  ans <- list(
    Sigmahat = sigma_hat, lambdahat = lambda_hat,
    Sigmasample = sample_covariance_matrix, Target = target,
    centered = centered
  )
  class(ans) <- "shrinkcovmathat"
  ans
}


shrinkcovmat_equal <- function(data, centered, p, n) {
  if (!centered) {
    if (n < 4) stop("The number of columns should be greater than 3")
    sample_covariance_matrix <- cov(t(data))
    trace_statistics <- trace_stats_uncentered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    nu_hat <- trace_sigma_hat / p
    trace_sigma_squared_hat <- trace_statistics[2]
    lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
      (n * trace_sigma_squared_hat + (p - n + 1) / p * trace_sigma_hat^2)
    lambda_hat <- min(lambda_hat, 1)
    } else {
      if (n < 2) stop("The number of columns should be greater than 1")
    sample_covariance_matrix <- tcrossprod(data) / n
    trace_statistics <- trace_stats_centered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    nu_hat <- trace_sigma_hat / p
    trace_sigma_squared_hat <- trace_statistics[2]
    lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
      ((n + 1) * trace_sigma_squared_hat + (p - n) / p * trace_sigma_hat^2)
    lambda_hat <- min(lambda_hat, 1)
  }
  if (lambda_hat < 1) {
    sigmahat <- (1 - lambda_hat) * sample_covariance_matrix +
      diag(nu_hat * lambda_hat, p)
  } else {
    sigmahat <- diag(lambda_hat * nu_hat, p)
  }
  target <- diag(nu_hat, p)
  ans <- list(
    Sigmahat = sigmahat, lambdahat = lambda_hat,
    Sigmasample = sample_covariance_matrix, Target = target,
    centered = centered
  )
  class(ans) <- "shrinkcovmathat"
  ans
}


shrinkcovmat_identity <- function(data, centered, p, n) {
  if (!centered) {  
  if (n < 4) stop("The number of columns should be greater than 3")
  sample_covariance_matrix <- cov(t(data))
  trace_statistics <- trace_stats_uncentered(data) # nolintr
  trace_sigma_hat <- trace_statistics[1]
  trace_sigma_squared_hat <- trace_statistics[2]
  lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
    (n * trace_sigma_squared_hat + trace_sigma_hat^2 -
       2 * trace_sigma_hat * (n - 1) + p * (n - 1))
  lambda_hat <- max(0, min(lambda_hat, 1))
  } else {
    if (n < 2) stop("The number of columns should be greater than 1")
    sample_covariance_matrix <- tcrossprod(data) / n
    trace_statistics <- trace_stats_centered(data) # nolintr
    trace_sigma_hat <- trace_statistics[1]
    trace_sigma_squared_hat <- trace_statistics[2]
    lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
      ((n + 1) * trace_sigma_squared_hat + trace_sigma_hat^2 -
         2 * trace_sigma_hat * n + p * n)
    lambda_hat <- max(0, min(lambda_hat, 1))
    }
if (lambda_hat < 1) {
  sigma_hat <- (1 - lambda_hat) * sample_covariance_matrix +
      diag(lambda_hat, p)
  } else {
    sigma_hat <- diag(lambda_hat, p)
  }
  target <- diag(p)
  ans <- list(
    Sigmahat = sigma_hat, lambdahat = lambda_hat,
    Sigmasample = sample_covariance_matrix, Target = target,
    centered = centered
  )
  class(ans) <- "shrinkcovmathat"
  ans
}
