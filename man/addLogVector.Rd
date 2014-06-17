\name{addLogVector}
\alias{addLogVector}
\title{
  Add a vector of log values.
}
\description{
  When a summation must be taken on a vector of numbers that are on the log scale, transforming them back to the non-log scale and taking the sum is slow.  This function takes the exp() of the values only when neccessary.
}
\usage{
  addLogVector(x)
}
\arguments{
  \item{x}{
    Numeric vector containing values on a log scale to be summed on the untransformed scale.
  }
}
\details{
  This function checks to see if the difference between the maximum values and the remaining values is less than the machine precision. If it is, then the exp() is taken for those values that differ by less than the machine precision, they are summed and returned to a log scale. If the maximum value is differs from the other values by greater than the machine precision, then return the maximum value.
}
\value{
  Numeric value containing the sum of the values on a log scale.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{addLog}}
}
\examples{
  addLogVector(log(1:10))
}
