shrinkcovmathat.default <- function(x, ...) {
    object <- list()
    object$Sigmahat <- x$Sigmahat
    object$lambdahat <- x$lambdahat
    object$Sigmasample <- x$Sigmasample
    object$centered <- x$centered
    object$Target <- x$Target
    class(object) <- "shrinkcovmathat"
    object
}
70
70
70
