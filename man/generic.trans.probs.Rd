\name{generic.trans.probs}
\alias{generic.trans.probs}
\title{
  Generic transition probabilities
}
\description{
  Generic function to provide starting transition probabilities. Not ready for use in this version.
}
\usage{
  generic.trans.probs(states, snps, chr = c(1:19, "X"), sex = c("M", "F"))
}
\arguments{
  \item{states}{
    Character vector containing the genotype states coded as letters.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
  \item{chr}{
    Character containing the chromosome to use.
  }
  \item{sex}{
    Character that is either FALSE or M indicating the sex. Only used on x chromosome.
  }
}
\details{
  This function is a placeholder for a more sophisticated algorithm to be developed later. It is an effort to provide a set of starting transition probabilities for mapping populations other than the DO.
}
\value{
  A matrix of transition probabilities between genotype states.
}
\note{
  This function is still in development.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{do.trans.probs}}, \code{\link{cc.trans.probs}}
}
\examples{
  \dontrun{ generic.trans.probs(states, snps, chr = 1, sex = "M") }
}
\keyword{ transition }
