\name{colSumsLog}
\alias{colSumsLog}
\alias{rowSumsLog}
\title{
  Sum columns of log transformed data.
}
\description{
  Given a matrix of log transformed values, sum the rows or columns on the untransformed scale.
}
\usage{
colSumsLog(logmat)
}
\arguments{
  \item{logmat}{
     Numeric matrix of natural log transformed values.
  }
}
\details{
  See \code{\link{addLog}}.
}
\value{
  Numeric vector with values summed on an untransformed scale.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{addLog}}
}
\examples{
  colSumsLog(matrix(log(runif(100)), nrow = 10, ncol = 10))
}
\keyword{ log }
