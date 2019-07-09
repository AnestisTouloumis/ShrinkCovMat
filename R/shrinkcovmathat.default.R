shrinkcovmathat.default <- function(x, ...) { #  nolint
  object <- list()
  object$Sigmahat <- x$Sigmahat # nolint
  object$lambdahat <- x$lambdahat
  object$Sigmasample <- x$Sigmasample # nolint
  object$centered <- x$centered
  object$Target <- x$Target # nolint
  class(object) <- "shrinkcovmathat"
  object
}
