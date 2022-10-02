#' @export

print.targetsel <- function(x, ...) {
  cat("ESTIMATED SHRINKAGE INTENSITIES WITH TARGET MATRIX THE", "\n")
  cat("Spherical matrix :", round(x$lambda_hat_spherical, 4), "\n")
  cat("Identity  matrix :", round(x$lambda_hat_identity, 4), "\n")
  cat("Diagonal  matrix :", round(x$lambda_hat_diagonal, 4), "\n")
  cat("\nSAMPLE VARIANCES", "\n")
  cat("Range   :", round(x$range, 4), "\n")
  cat("Average :", round(x$average, 4), "\n")
}
