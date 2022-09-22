shrinkcovmat_unequal <- function(data, centered, p, n) {
  sample_covariance_matrix <- calculate_covariance_matrix(data, centered, n)
  trace_statistics <- calculate_trace_statistic(data, centered)
  trace_sigma_hat <- trace_statistics[1]
  trace_sigma_squared_hat <- trace_statistics[2]
  trace_diagonal_sigma_sq_hat <- trace_statistics[3]
  sample_variances <- calculate_sample_variances(data, centered)
  sample_size <- ifelse(centered, n + 1, n)
  lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat -
    (2 - 2 / sample_size) * trace_diagonal_sigma_sq_hat) /
    (sample_size * trace_sigma_squared_hat + trace_sigma_hat^2 -
      (sample_size + 1 - 2 / sample_size) * trace_diagonal_sigma_sq_hat)
  lambda_hat <- max(0, min(lambda_hat, 1))
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
  sample_covariance_matrix <- calculate_covariance_matrix(data, centered, n)
  trace_statistics <- calculate_trace_statistic(data, centered)
  trace_sigma_hat <- trace_statistics[1]
  trace_sigma_squared_hat <- trace_statistics[2]
  sample_size <- ifelse(centered, n + 1, n)
  nu_hat <- trace_sigma_hat / p
  lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
    (sample_size * trace_sigma_squared_hat + (p - sample_size + 1) / p * trace_sigma_hat^2)
  lambda_hat <- min(lambda_hat, 1)
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
  sample_covariance_matrix <- calculate_covariance_matrix(data, centered, n)
  trace_statistics <- calculate_trace_statistic(data, centered)
  trace_sigma_hat <- trace_statistics[1]
  trace_sigma_squared_hat <- trace_statistics[2]
  sample_size <- ifelse(centered, n + 1, n)
  lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
    (sample_size * trace_sigma_squared_hat + trace_sigma_hat^2 -
      2 * trace_sigma_hat * (sample_size - 1) + p * (sample_size - 1))
  lambda_hat <- max(0, min(lambda_hat, 1))
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


calculate_covariance_matrix <- function(data, centered, n) {
  if (centered) {
    ans <- tcrossprod(data) / n
  } else {
    ans <- cov(t(data))
  }
  ans
}


calculate_trace_statistic <- function(data, centered) {
  if (centered) {
    ans <- trace_stats_centered(data)
  } else {
    ans <- trace_stats_uncentered(data)
  }
  ans
}


calculate_sample_variances <- function(data, centered) {
  if (centered) {
    ans <- apply(data, 1, function(x) mean(x^2))
  } else {
    ans <- apply(data, 1, var)
  }
  ans
}
