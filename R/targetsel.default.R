targetsel.default <- function(x, ...) {
  object <- list()
  object$lambda_hat_spherical <- x$lambda_hat_spherical
  object$lambda_hat_identity <- x$lambda_hat_identity
  object$lambda_hat_diagonal <- x$lambda_hat_diagonal
  object$range <- x$range
  object$average <- x$average
  class(object) <- "targetsel"
  object
}
