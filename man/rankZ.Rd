\name{rankZ}
\alias{rankZ}
\title{
  Rank Z transformation
}
\description{
  This is also called the inverse normal transformation. We rank the data in x, divide by n-1, and take quantiles from the normal distribution using qnorm.
}
\usage{
  rankZ(x)
}
\arguments{
  \item{x}{
    Numeric vector of values to be transformed.
  }
}
\details{
  We often use this when there are hundreds or thousands of phenotypes. This allows us to calculate permutation derived thresholds once and use them for all phenotypes.
}
\value{
  Numeric vector with a normal distribution.
}
\author{
  Daniel Gatti
}
\examples{
  rankZ(rlnorm(20))
}
\keyword{ rankZ }
