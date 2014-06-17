\name{convert.allele.calls}
\alias{convert.allele.calls}
\title{
  Convert allele calls to numeric values.
}
\description{
  Converts allele calls in A,C,G,T,H,N format into numbers with 0: homozygous A, 1: heterozygous, 2: homozygous B, 3: no call.
}
\usage{
  convert.allele.calls(geno)
}
\arguments{
  \item{geno}{
    Character matrix containing A, C, G, T, H or N.
  }
}
\value{
  Numeric matrix containing 0, 1, 2 or 3.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ convert.allele.calls(geno) }
}
\keyword{ alleles }
