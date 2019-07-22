targetsel.default <- function(x, ...) { #  nolint
  object <- list()
  object$optimal_sphericity <- x$optimal_sphericity
  object$optimal_identity <- x$optimal_identity
  object$optimal_diagonal <- x$optimal_diagonal
  object$range <- x$range
  object$average <- x$average
  class(object) <- "targetsel"
  object
}
