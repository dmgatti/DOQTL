\name{summarize.genotype.transitions}
\alias{summarize.genotype.transitions}
\alias{summarize.by.samples}
\alias{summarize.by.snps}
\alias{num.recomb.plot}
\title{Summarize the genotype data output by the genotyping HMM.}
\description{
  These functions read in all of the individual genotype data files and summarizes the founder allele frequency by sample and SNP.
}
\usage{
  summarize.genotype.transitions(path = ".", snps)
  summarize.by.snps(path = ".", snps)
  summarize.by.samples(path = ".", snps)
  num.recomb.plot(results, gen) 
}
\arguments{
  \item{path}{
    Character, full path to the genotype directory where the *.genotype.probs.Rdata files are stored.
  }
  \item{snps}{
     Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
  \item{results}{
    Data.frame containing results as summarized by summarize.genotype.transitions.
  }
  \item{gen}{
    Numeric vector containing DO outbreeding generation for all samples in results. Must have sample IDs in names.
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
\keyword{ Diversity Outbred }
\keyword{ genotyping }
\keyword{ HMM }
