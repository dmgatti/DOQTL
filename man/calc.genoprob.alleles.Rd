\name{calc.genoprob.alleles}
\alias{calc.genoprob.alleles}
\title{Calculate the founder genotype probabilities at each SNP using allele calls.}
\description{
  This function performs genome reconstruction using allele calls. We recommend using allele intensities where available because they often produce better genotype reconstructions.
}
\usage{
  calc.genoprob.alleles(data, chr, founders, snps, output.dir = ".", 
  trans.prob.fxn = do.trans.probs, plot = FALSE)
}
\arguments{
  \item{data}{
    A list with named elements containing the information needed to reconstruct genomes. 
	
       geno: Character matrix, num.samples x num.snps, with allele calls (A,C,G,T,H or N) for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   gen: Character matrix containing the generation of DO outbreeding for each sample. For the DO, this should be "DO" followed by a number with no space between them. For CC mice, this should be CC. Sample IDs must be in names.	   
  }
  \item{chr}{
     Character vector containing chromosomes to run. Must match the chromosome IDs in the snps table.  "all" (default) will run all chromosomes.
  }
  \item{founders}{
     List containing founder information for non-DO or CC crosses.
    When method = allele:
       geno: Character matrix, num.samples x num.snps, with allele calls (A,C,G,T,H or N) for all founders and FALSE1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   code: Character vector containing two letter genotype codes for each founder sample. Sample IDs must be in names.
  }
  \item{snps}{
     Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively. \emph{Not required for DO.}
  }
  \item{output.dir}{
     Character string containing the full path where output should be written. The directory must exist already.
  }
  \item{plot}{
     Boolean that is true if the user would like to plot a sample chromosome as the model progresses. Default  =  TRUE.
  }
  \item{trans.prob.fxn}{
     Function to use when computing the transtion probabilities between markers.
  }
}
\value{
  No value is returned.  The output files are written to output.dir.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{
    calc.genoprob.alleles(data, chr = 1:19, founders = founders, 
	snps = snps, output.dir = "do.data")
  }
}
\keyword{ MUGA }
\keyword{ genotyping }
\keyword{ HMM }
