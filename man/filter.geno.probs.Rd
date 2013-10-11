\name{filter.geno.probs}
\alias{filter.geno.probs}
\title{Remove SNPs where the genotype probabilities are too low for one founder state}
\description{
  This function accetps the 8 state founder probabilities and searches for SNPs where the founder probabilities are too low for one founder.  It removes these SNPs.  The QTL mapping model becomes numerically unstable if these SNPs are not removed.
}
\usage{
  filter.geno.probs(geno)
}
\arguments{
  \item{geno}{
    Numeric 3D array, with samples in dim[1], states = dim[2] and SNPs in dim[3].
  }
}
\value{
Numeric 3D array, with samples in dim[1], states = dim[2] and SNPs in dim[3].  The SNPs with low probabilities for a single founder have been removed.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{scanone}}, \code{\link{scanone.perms}}
}
\examples{
  \dontrun{
    filter.geno.probs(geno)
  }
}
\keyword{ QTL }
\keyword{ genotypes }

