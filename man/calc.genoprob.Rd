\name{calc.genoprob}
\alias{calc.genoprob}
\title{Calculate the founder genotype probabilities at each SNP.}
\description{
  This function performs genome reconstruction using either allele calls or allele intensities. We recommend using allele intensities where available because they often produce better genotype reconstructions.
}
\usage{
  calc.genoprob(data, chr = "all", output.dir = ".", plot = TRUE, 
  array = c("gigamuga", "megamuga", "muga", "other"), 
  sampletype = c("DO", "CC", "DOF1", "other"), method = c("intensity", "allele"), 
  founders, transprobs, snps)
}
\arguments{
  \item{data}{
    A list with named elements containing the information needed to reconstruct genomes. 
	
    When method = intensity:
       x: Numeric matrix, num.samples x num.snps, with X intensities for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or "F" indicating sex. Sample IDs must be in names.
	   gen: Character matrix containing the generation of DO outbreeding for each sample. For the DO, this should be "DO" followed by a number with no space between them. For CC mice, this should be CC. Sample IDs must be in names.
	   
    When method = allele:
       geno: Character matrix, num.samples x num.snps, with allele calls (A,C,G,T,H or N) for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   gen: Character matrix containing the generation of DO outbreeding for each sample. For the DO, this should be "DO" followed by a number with no space between them. For CC mice, this should be CC. Sample IDs must be in names.	   
  }
  \item{chr}{
     Character vector containing chromosomes to run. Must match the chromosome IDs in the snps table.  "all" (default) will run all chromosomes.
  }
  \item{output.dir}{
     Character string containing the full path where output should be written. The directory must exist already.
  }
  \item{plot}{
     Boolean that is true if the user would like to plot a sample chromosome as the model progresses. Default  =  TRUE.
  }
  \item{array}{
     Character string indicating the array type. Must be one of "gigamuga", "megamuga", "muga" or "other". Default equals "gigamuga".
  }
  \item{sampletype}{
     Character string containing the type of samples being run. Must be one of "DO", "CC", "DOF1", or "other". Default equals "DO".
  }
  \item{method}{
     Character string containing method of genome reconstruction. Must be one of "intensity" or "allele". Default equals "intensity".
  }
  \item{founders}{
     List containing founder information for non-DO or CC crosses. \emph{Not required for DO.}
	 
    When method = intensity:
       x: Numeric matrix, num.samples x num.snps, with X intensities for all founders and F1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all founders and F1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or "F" indicating sex. Sample IDs must be in names.
	   code: Character vector containing two letter genotype codes for each founder sample. Sample IDs must be in names.
	   
    When method = allele:
       geno: Character matrix, num.samples x num.snps, with allele calls for all founders and F1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or "F" indicating sex. Sample IDs must be in names.
	   code: Character vector containing two letter genotype codes for each founder sample. Sample IDs must be in names.
	   
    When sampletype = DOF1
	   x: Numeric matrix, num.samples x num.snps, with X intensities for all founders and F1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all founders and F1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or "F" indicating sex. Sample IDs must be in names.
	   code: Character vector containing two letter genotype codes for each founder sample. This should be "II" for an inbred mutant strain. Sample IDs must be in names.
	   direction: Character string that is either "DOxMUT" if a female DO was crossed with a male mutant mouse or "MUTxDO" if a female mutant mouse was crossed with a male DO. This affects the genotyping of teh X chromosome.
  }
  \item{transprobs}{
     Function to call to estimate the transition probabilities between markers for non-DO samples. \emph{Not required for DO.}
  }
  \item{snps}{
     Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively. \emph{Not required for DO.}
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
    calc.genoprob(cross, chr.to.run = 1:19, output.dir = "do.data", plot = FALSE, 
	init.means = NULL, init.covars = NULL)
  }
}
\keyword{ MUGA }
\keyword{ genotyping }
\keyword{ HMM }
