\name{create.genotype.states}
\alias{create.genotype.states}
\title{
  Create genotype states.
}
\description{
  Given a set of founders, create all of the possible unphased genotype states between them.
}
\usage{
  create.genotype.states(founders, sampletype = c("DO", "CC", "HS", "HSrat"))
}
\arguments{
  \item{founders}{
     Character vector of letter codes indicating the founders.
  }
  \item{sampletype}{
     Character string that is one of DO, HS, CC or HSrat.
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
  create.genotype.states(founders = LETTERS[1:8], sampletype = "DO")
}
\keyword{ genotype }
\keyword{ HMM }
