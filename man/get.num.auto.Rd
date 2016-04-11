\name{get.num.auto}
\alias{get.num.auto}
\title{
  Get the number of autosomes
}
\description{
  Given a data.frame of SNPs, return the number of autosomes by counting the number of numeric chromosomes.
}
\usage{
  get.num.auto(snps)
}
\arguments{
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
}
\value{
  Numeric value containing the number of autosomes.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ get.num.auto(snps) }
}
\keyword{ autosomes }
