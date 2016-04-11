\name{genome.summary.plots}
\alias{genotype.by.sample.barplot}
\alias{genotype.by.snp.barplot}
\title{
  Genome summary plots
}
\description{
  These plots show summaries of the founder haplotype proportions across SNPs or samples.
}
\usage{
  genotype.by.sample.barplot(results)
  genotype.by.snp.barplot(results, snps)
}
\arguments{
  \item{results}{
    Data.frame as produced by \code{\link{summarize.by.snps}} or \code{\link{summarize.by.samples}}.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
}
\details{
  Barplots across samples or SNPs are produced from the relevant summary files. This is intended to be run after DO genomes have been reconstructed for a large set of samples to verify that founder allele frequencies are consistent across SNPs and samples.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{summarize.by.snps}}, \code{\link{summarize.by.samples}}
}
\examples{
  \dontrun{ 
    load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
    results = summarize.by.samples(path = ".", snps = muga_snps)
    genotype.by.sample.barplot(results) 
  }
}
\keyword{ summary }
\keyword{ genome }
