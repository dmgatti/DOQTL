\name{get.trans.probs}
\alias{get.trans.probs}
\title{
  Get the transition probabilities between markers.
}
\description{
  Based on the DO founders and the recombination fraction, calculate the transition probability between markers on one chromosome.
}
\usage{
get.trans.probs(r, do.gen, alpha, chr = c(1:19, "X"), sex = c("M", "F"))
}
\arguments{
  \item{r}{
    Double vector of recombination fractions between SNPs.
  }
  \item{do.gen}{
    Integer containing the outbreeding generation of DO.
  }
  \item{alpha}{
    Double vector, proportion of preCC progenitors at generation k. Generation numbers in the names.
  }
  \item{chr}{
    Character, one of 1:19, X.
  }
  \item{sex}{
    Character, either M or FALSE. Only used on X chromosome.
  }
}
\value{
   Numeric three dimensional array of transition probabilities between each pair of markers (num.states by num.states by num.markers - 1).
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
  Daniel Gatti
}
\examples{
  \dontrun{ get.trans.probs(r, do.gen, alpha, chr = c(1:19, "X"), sex = c("M", "F")) }
}
\keyword{ transition }
\keyword{ HMM }
