\name{get.max.geno}
\alias{get.max.geno}
\title{
  Get the genotype with the highest probability
}
\description{
  For each sample at each marker, return the genotype state with the highest probability.
}
\usage{
  get.max.geno(probs)
}
\arguments{
  \item{probs}{
    Three dimensional numeric matrix containing the founder haplotype contributions.
  }
}
\details{
  If the maximum probability is greater than 0.75, it is called homozygous. Otherwise, we take the two highest probabilities and call heterozygous.
}
\value{
  Character matrix containing the genotype with the highest probability at each locus.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ get.max.geno(probs) }
}
\keyword{ genotype }
