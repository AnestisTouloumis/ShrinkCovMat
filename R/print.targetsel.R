#' @export

print.targetsel <- function(x, ...) {
  cat("OPTIMAL SHRINKAGE INTENSITIES FOR THE TARGET MATRIX WITH", "\n")
  cat("Equal variances   :", round(x$optimal_sphericity, 4), "\n")
  cat("Unit variances    :", round(x$optimal_identity, 4), "\n")
  cat("Unequal variances :", round(x$optimal_diagonal, 4), "\n")
  cat("\nSAMPLE VARIANCES", "\n")
  cat("Range   :", round(x$range, 4), "\n")
  cat("Average :", round(x$average, 4), "\n")
}
