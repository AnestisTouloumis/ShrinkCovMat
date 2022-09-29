calculate_sample_covariance_matrix <- function(data, centered, n) {
  if (!centered) {
    ans <- cov(t(data))
  } else {
    ans <- tcrossprod(data) / n
  }
  ans
}


calculate_trace_statistics <- function(data, centered) {
  if (!centered) {
    ans <- trace_stats_uncentered(data)
  } else {
    ans <- trace_stats_centered(data)
  }
  ans
}


calculate_sample_variances <- function(data, centered) {
  if (!centered) {
    ans <- apply(data, 1, var)
  } else {
    ans <- apply(data, 1, function(x) mean(x^2))
  }
  ans
}


calculate_lambda_hat <- function(trace_statistics, sample_size, p, target) {
  trace_sigma_hat <- trace_statistics[1]
  trace_sigma_squared_hat <- trace_statistics[2]
  if (target == "diagonal") {
    trace_diagonal_sigma_sq_hat <- trace_statistics[3]
    ans <- (trace_sigma_hat^2 + trace_sigma_squared_hat -
      (2 - 2 / sample_size) * trace_diagonal_sigma_sq_hat) /
      (sample_size * trace_sigma_squared_hat + trace_sigma_hat^2 -
        (sample_size + 1 - 2 / sample_size) * trace_diagonal_sigma_sq_hat)
  } else if (target == "spherical") {
    ans <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
      (sample_size * trace_sigma_squared_hat + (p - sample_size + 1) / p * trace_sigma_hat^2)
  } else {
    ans <- (trace_sigma_hat^2 + trace_sigma_squared_hat) /
      (sample_size * trace_sigma_squared_hat + trace_sigma_hat^2 -
        2 * trace_sigma_hat * (sample_size - 1) + p * (sample_size - 1))
  }
  ans <- max(0, min(ans, 1))
  ans
}


calculate_target_matrix <- function(data, centered, p, target) {
  if (target == "diagonal") {
    sample_variances <- calculate_sample_variances(data, centered)
    ans <- diag(sample_variances, p)
  } else if (target == "spherical") {
    sample_variances <- calculate_sample_variances(data, centered)
    ans <- diag(mean(sample_variances), p)
  } else {
    ans <- diag(p)
  }
  ans
}


calculate_shrinkage_covariance_matrix <- function(sample_covariance_matrix,
                                                  lambda_hat, target_diagonal) {
  if (lambda_hat < 1) {
    ans <- (1 - lambda_hat) * sample_covariance_matrix
    diag(ans) <- diag(ans) + lambda_hat * target_diagonal
  } else {
    ans <- diag(lambda_hat * target_diagonal)
  }
  ans
}
