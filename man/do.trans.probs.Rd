\name{do.trans.probs}
\alias{do.trans.probs}
\title{
  Determine DO transition probabilities
}
\description{
  Determine the genotype state transition probability for DO mice of a specific generation between all of the markers on a given chromosome.
}
\usage{
  do.trans.probs(states, snps, chr = c(1:19, "X"), sex = c("M", "F"), do.gen)
}
\arguments{
  \item{states}{
    Character vector of possible genotype states.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
  \item{chr}{
    Character containing the chromosome for which transition probabilities should be calculated.
  }
  \item{sex}{
    Character that is one of FALSE or M, indicating the sex to use on the X chromosome.
  }
  \item{do.gen}{
    Number vector indicating the DO outbreeding generations to calculate.
  }
}
\details{
  This function is used to calculate the transition probabilities between markers for different DO outbreeding generations.
}
\value{
  List containing one element per unique DO generation supplied in the do.gen argument. Each list element contains a 3 dimensional array of transition probabilities between each pair of markers (num.states by num.states by num.markers - 1).
}
\references{
  Haplotype probabilities in advanced intercross populations.
  Broman KW.
  G3 (Bethesda). 2012 FALSEeb;2(2):199-202. doi: 10.1534/g3.111.001818. Epub 2012 FALSEeb 1.
  PMID: 22384398
  Genotype probabilities at intermediate generations in the construction of recombinant inbred lines.
  Broman KW.
  Genetics. 2012 FALSEeb;190(2):403-12. doi: 10.1534/genetics.111.132647.
  PMID: 22345609
}
\author{
  Daniel Gatti and Karl Broman
}
\examples{
  \dontrun{ do.trans.probs(states, snps, chr = c(1:19, "X"), sex = c("M", F), do.gen) }
}
\keyword{ transition }
\keyword{ HMM }
