\name{condense.model.probs}
\alias{condense.model.probs}
\alias{get.additive}
\alias{get.dominance}
\alias{get.full}
\title{
  Condense 36 state genotypes down to founder genotypes.
}
\description{
  Additive condenses the heterozygous genotype calls down to the founder allele contributions. Dominance will eventually provide additive and dominance values (method currently uncertain). FALSEull simply gathers all of the genotype probabilities together. 
}
\usage{
  condense.model.probs(path = ".", write, model = c("additive", "dominance", "full"), cross = "DO")
  get.additive(files, samples)
  get.dominance(files, samples)
  get.full(files, samples)
}
\arguments{
  \item{path}{
    Character containing the path to the *.Rdata files that contain the genotype probabilities.
  }
  \item{write}{
    Character that is the filename to write the results to.
  }
  \item{model}{
    Character string that is one of "additive", "dominance" or "full". See details.
  }
  \item{cross}{
    Character string containing the type of cross. Typically "DO, "CC", "HS" or "DOF1". This is used in downsteam analyses to produce the coefficient plots and select strains for association mapping.
  }
  \item{files}{
    Vector of files to read.
  }
  \item{samples}{
    Character vector of sample IDs that match the files argument.
  }
}
\details{
  get.additive, get.dominance and get.full are helper functions.
}
\value{
  Three dimensional array of haplotype or genotype probabilities. Num. samples by num. founders (or genotypes) by num. SNPs. Writes out to a *.Rdata file.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{calc.genoprob}}
}
\examples{
  \dontrun{ condense.model.probs(write = "model.probs.Rdata") }
}
\keyword{ HMM }
