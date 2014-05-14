\name{addLog}
\alias{addLog}
\title{
  Add two log values.
}
\description{
  When two numbers that are on the log scale must be summed, transforming them back to the non-log scale and taking the sum is slow.  This function takes the exp() of the values only when neccessary.
}
\usage{
  addLog(x, y)
}
\arguments{
  \item{x}{
    Numeric containing a value on the log scale.
  }
  \item{y}{
    Numeric containing a value on the log scale.
  }
}
\details{
  This function checks to see if the difference between the maximum value and the other value is less than the machine precision. If it is, then the exp() is taken for those values that differ by less than the machine precision, they are summed and returned to a log scale. If the maximum value is differs from the other values by greater than the machine precision, then return the maximum value.
}
\value{
  Numeric value containing the sum of the arguments on a log scale.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{addLogVector}}
}
\examples{
  addLog(log(10), log(1))
}
