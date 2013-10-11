\name{argmax.geno}
\alias{argmax.geno}
\title{Find the most probable genotype for all samples using the Viterbi algorithm.}
\description{
  This function reads in the X and Y intensity values for all samples as well as the final means and covariances for the clusters at each SNP.  It then constructs emission and transtion probability matrices at each SNP and uses the Viterbi algorithm to find the more probable genotype.
}
\usage{
  argmax.geno(X.file, Y.file, sex.gen.file, mean.files, covar.files, chr, snps, direction = c("forward", "backward"))
}
\arguments{
  \item{X.file}{
    character, required: the full path to the X intensity file for all samples.
  }
  \item{Y.file}{
    character, required: the full path to the Y intensity file for all samples.
  }
  \item{sex.gen.file}{
    character, required: the full path to the a file containing the sex and DO outbreeding generation, with sample IDs in rownames.
  }
  \item{mean.files}{
    character vector, required: the full path to the mean files output by calc.genoprob. These will typically be *.Rdata files with one file per chromosome.
  }
  \item{covar.files}{
    character vector, required: the full path to the covariance files output by calc.genoprob. These will typically be *.Rdata files with one file per chromosome.
  }
  \item{chr}{
    character vector, required: the chromosomes to process.
  }
  \item{snps}{
    data.frame, required: SNPs to be used in 4 columns, with SNP ID, Chr, Mb and cM, in that order.
  }
  \item{direction}{
    character, required: one of "forward" or "backward".  Forward denotes that the algorithm will start at the proximal end of the chromosome. Backward indicate that it will start at the distal end.
  }
}
\value{
  List with one element per chromosome.  Each element is a character matrix with two letter, unphased genotype calls.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{calc.genoprob}}, \code{\link{plot.genotype}}
}
\examples{
  \dontrun{max.geno.back = argmax.geno(X.file = "/demo/MUGA/x.filt.batch.norm.txt",
           Y.file = "/demo/MUGA/y.filt.batch.norm.txt",
           sex.gen.file = "/demo/MUGA/demo.sex.gen.txt",
           mean.files = dir(path = "/demo/HMM", pattern = "means.Rdata", full.names = T),
           covar.files = dir(path = "/demo/HMM", pattern = "covars.Rdata", full.names = T),
           chr = c(1:19, "X"), snps, direction = "forward")}
}
\keyword{ Viterbi }
\keyword{ genotype }
