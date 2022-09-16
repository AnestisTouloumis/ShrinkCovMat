#' Linear Shrinkage of the Sample Covariance
#'
#' Provides a nonparametric Stein-type shrinkage estimator of the covariance
#' matrix.
#'
#' The rows of the data matrix \code{data} correspond to variables and the
#' columns to subjects.
#'
#' @param data a numeric matrix containing the data.
#' @param centered a logical indicating if the mean vector is the zero vector.
#' @param target a character indicating the target matrix. Options include 
#' 'sphericity', 'identity' or 'diagonality'. 
#' @return Returns an object of the class 'shrinkcovmathat' that has
#' components: \item{Sigmahat}{The Stein-type shrinkage estimator of the
#' covariance matrix.} \item{lambdahat}{The estimated optimal shrinkage
#' intensity.} \item{Sigmasample}{The sample covariance matrix.}
#' \item{Target}{The target covariance matrix.} \item{centered}{If the data are
#' centered around their mean vector.}
#' @author Anestis Touloumis
#' @seealso \code{\link{targetselection}}.
#' @references Touloumis, A. (2015) nonparametric Stein-type Shrinkage
#' Covariance Matrix Estimators in High-Dimensional Settings.
#' \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
#' @examples
#' data(colon)
#' normal_group <- colon[, 1:40]
#' tumor_group <- colon[, 41:62]
#' targetselection(normal_group)
#' sigma_hat_normal_group <- shrinkcovmat(normal_group)
#' sigma_hat_normal_group
#' targetselection(normal_group)
#' sigma_hat_tumor_group <- shrinkcovmat(tumor_group)
#' sigma_hat_tumor_group
#' @export
shrinkcovmat <- function(data, centered = FALSE, target = "sphericity") { # nolint
  targets <- c("sphericity", "identity", "diagonality")
  if (!is.element(target, targets)) {
    stop("'target' must be 'sphericity', 'identity' or 'diagonality'")
    }
  if (target == "sphericity") {
    ans <- shrinkcovmat.equal(data, centered)
    } else if (target == "diagonality") {
      ans <- shrinkcovmat.unequal(data, centered)
      } else {
        ans <- shrinkcovmat.identity(data, centered)
        }
  ans
}