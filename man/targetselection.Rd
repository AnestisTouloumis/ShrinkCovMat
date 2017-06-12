\name{targetselection}
\alias{targetselection}
\title{
Target Matrix Selection
}
\description{
Implements the rule of thumb proposed by Touloumis (2015) for target matrix selection. If the estimated optimal shrinkage intensities of the three target matrices are of similar magnitude, then the average and the range of the sample variances should be inspected in order to adopt the most plausible target matrix.
}
\usage{
targetselection(data, centered = FALSE)
}
\arguments{
  \item{data}{
a numeric matrix containing the data.
}
  \item{centered}{
a logical indicating if the mean vector is the zero vector.
}
}
\details{
The rows of the data matrix \code{data} correspond to variables and the columns to subjects. 
}
\value{
Prints the estimated optimal shrinkage intensities and the range and the average of the sample variances.
}
\references{
Touloumis, A. (2015) Nonparametric Stein-type Shrinkage Covariance Matrix Estimators in High-Dimensional Settings. \emph{Computational Statistics & Data Analysis} \bold{83}, 251--261. 
}
\author{
Anestis Touloumis
}
\examples{
data(colon)
normal.group <- colon[,1:40]
targetselection(normal.group)
## Similar intensities, the range of the sample variances is small 
## and the average is not close to one. The scaled identity matrix 
## seems to be the most suitable target matrix for the normal group

colon.group <- colon[,41:62]
targetselection(colon.group)
## Similar intensities, the range of the sample variances is small 
## and the average is not close to one. The scaled identity matrix 
## seems to be the most suitable target matrix for the colon group
}