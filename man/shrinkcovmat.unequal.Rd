\name{shrinkcovmat.unequal}
\alias{shrinkcovmat.unequal}
\title{
Shrinking the Sample Covariance Matrix Towards a Diagonal Matrix with Diagonal Elements the Sample Variances.
}
\description{
Provides a nonparametric Stein-type shrinkage estimator of the covariance matrix that is a linear combination of the sample covariance matrix and of the diagonal matrix with elements the corresponding sample variances on the diagonal and zeros elsewhere.
}
\usage{
shrinkcovmat.unequal(data, centered = FALSE)
}
\arguments{
  \item{data}{
a numeric matrix containing the data.
}
  \item{centered}{
a logical indicating if the vectors are centered around their mean vector.
}
}
\details{
The rows of the data matrix \code{data} correspond to variables and the columns to subjects. 
}
\value{
Returns an object of the class "shrinkcovmathat" that has components:
\item{Sigmahat}{The Stein-type shrinkage estimator of the covariance matrix.}
\item{lambdahat}{The estimated optimal shrinkage intensity.}
\item{Sigmasample}{The sample covariance matrix.}
\item{Target}{The target covariance matrix.}
\item{centered}{If the data are centered around their mean vector.}
}
\references{
Touloumis, A. (2015) Nonparametric Stein-type Shrinkage Covariance Matrix Estimators in High-Dimensional Settings. \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261.
}
\author{
Anestis Touloumis
}
\seealso{
\code{\link{shrinkcovmat.equal}} and \code{\link{shrinkcovmat.identity}}.
}
\examples{
data(colon)
normal.group <- colon[, 1:40]
colon.group <- colon[, 41:62]
Sigmahat.normal <- shrinkcovmat.unequal(normal.group)
Sigmahat.normal
Sigmahat.colon <- shrinkcovmat.unequal(colon.group)
Sigmahat.colon
}

