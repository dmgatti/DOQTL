\name{create.genotype.states}
\alias{create.genotype.states}
\title{
  Create genotype states.
}
\description{
  Given a set of founders, create all of the possible unphased genotype states between them.
}
\usage{
  create.genotype.states(founders)
}
\arguments{
  \item{founders}{
     Character vector of letter codes indicating the founders.
  }
}
\details{
  Given a set of founder IDs, create all possible unphased genotypes and sort them.
}
\value{
  Character vector of unphased genotypes that can be created from the given founders.
}
\author{
  Daniel Gatti
}
\examples{
  create.genotype.states(founders = LETTERS[1:8])
}
\keyword{ genotype }
\keyword{ HMM }
