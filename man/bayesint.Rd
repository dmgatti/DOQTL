\name{bayesint}
\alias{bayesint}
\title{Find a Bayesian Credible Interval around a QTL.}
\description{
  This function normalizes the area under the QTL curve on the given chromosome and finds a region that is 95% of the area.  Other probabilities can also be given.
}
\usage{
  bayesint(qtl, chr, prob = 0.95, expandtomarkers = TRUE)
}
\arguments{
  \item{qtl}{
    data.frame: four columns with SNP ID, Chr, position, and LOD in each column.
  }
  \item{chr}{
    character: the chromosome on which the QTL lies.
  }
  \item{prob}{
    numeric: must be between 0 and 1
  }
  \item{expandtomarkers}{
    boolean: if TRUE, expand the QTL interval to the nearest flanking markers.  Default = TRUE.
  }
}
\value{
  Data frame with the SNP ID, Chr, position and LOD for the left and right side of the interval and the maximum QTL.
}
\author{
  Daniel Gatti
}
\references{
  Saunak, S. (2001) A Statistical Framework for Quantitative Trait Mapping.  \emph{Genetics}, \bold{159} (1), 371--387.
}
\seealso{
  \code{\link{scanone}}, \code{\link{scanone.perm}}
}
\examples{
  \dontrun{bayesint(qtl, 1)}
}
\keyword{ QTL }
