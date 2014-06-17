\name{convert.genotypes}
\alias{convert.genotypes}
\title{
  Convert the genotype data from A,C,G,T format to A, H, B, N.
}
\description{
  Convert the genotype data from A,C,G,T format to A, H, B, N.
}
\usage{
  convert.genotypes(geno)
}
\arguments{
  \item{geno}{
     Character matrix containing A, C, G, T allele calls.
  }
}
\value{
  Character matrix containing A, H, B, N.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ convert.genotypes(geno) }
}
\keyword{ alleles }
