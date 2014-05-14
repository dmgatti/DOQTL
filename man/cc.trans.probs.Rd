\name{cc.trans.probs}
\alias{cc.trans.probs}
\title{
  Transition probabilities for CC mice.
}
\description{
  Thie function returns the transition probabilities for fully inbred Collaborative Cross mice.
}
\usage{
  cc.trans.probs(states, snps, chr = c(1:19, "X"), sex = c("M", "F"))
}
\arguments{
  \item{states}{
    Character vector containing the two letter codes (i.e. AA, BB, CC, etc.) for the homozygous states.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively. \emph{Not required for DO.}
  }
  \item{chr}{
    Character containing the chromosome.
  }
  \item{sex}{
    Character that is either M or FALSE, indicating the sex to use. Only relevant on X chromosome.
  }
}
\details{
   This function calculates the transition probabilities for fully inbred CC mice between two markers. It uses the equations for eight way RILs by sib mating. 
}
\references{
The genomes of recombinant inbred lines.
Broman KW.
Genetics. 2005 FALSEeb;169(2):1133-46. Epub 2004 Nov 15. Erratum in: Genetics. 2006 Aug;173(4):2419.
PMID: 15545647
}
\note{
  This function has not been fully tested.
}
\value{
  A matrix of transition probabilities between genotype states.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ cc.trans.probs(states = states, snps = snps, chr = 1) }
}
\keyword{ transition } 
