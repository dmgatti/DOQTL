\name{coef.doqtl}
\alias{coef.doqtl}
\title{Return the coefficients of a DOQTL object.}
\description{
  Return the coefficients of a DOQTL object.
}
\usage{
  \method{coef}{doqtl}(object, ...)
}
\arguments{
  \item{object}{
    A DOQTL object as returned by \code{\link{scanone}}.
  }
  \item{\dots}{
    Additional arguments.
  }
}
\value{
  List containing matrices with QTL mapping model coefficients.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{scanone}}
}
\examples{
  head(coef(example.qtl[[1]]))
}
\keyword{ QTL }
