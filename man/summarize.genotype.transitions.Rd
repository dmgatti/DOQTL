\name{summarize.genotype.transitions}
\alias{summarize.genotype.transitions}
\title{Summarize the genotype data output by the genotyping HMM.}
\description{
  This function reads in all of the individual genotype data files and summarizes the founder allele frequency by sample and SNP.
}
\usage{
  summarize.genotype.transitions(path = ".", snps)
}
\arguments{
  \item{path}{
    Character, full path to the genotype directory where the *.genotype.probs.Rdata files are stored.
  }
  \item{snps}{
     Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
}
\value{
  Data.frame with six columns: sample, proximal genotype, distal genotype, chr, proximal Mb, distal Mb.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ geno.sum = summarize.genotype.transitions() }
}
\keyword{ Diverisity Outbred }
\keyword{ genotyping }
\keyword{ HMM }

