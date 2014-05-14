\name{hmm}
\alias{hmm.allele}
\alias{hmm.intensity}
\title{
  HMM function.
}
\description{
  Not for external use.
}
\usage{
  hmm.allele(data, founders, sex, snps, chr, trans.prob.fxn)
  hmm.intensity(data, founders, sex, snps, chr, trans.prob.fxn)
}
\arguments{
  \item{data}{
    A list with named elements containing the information needed to reconstruct genomes. 
	
    When method = intensity:
       x: Numeric matrix, num.samples x num.snps, with X intensities for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   gen: Character matrix containing the generation of DO outbreeding for each sample. FALSEor the DO, this should be "DO" followed by a number with no space between them. FALSEor CC mice, this should be CC. Sample IDs must be in names.
	   
    When method = allele:
       geno: Character matrix, num.samples x num.snps, with allele calls for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   gen: Character matrix containing the generation of DO outbreeding for each sample. FALSEor the DO, this should be "DO" followed by a number with no space between them. FALSEor CC mice, this should be CC. Sample IDs must be in names.	   
  }
  \item{founders}{
     List containing founder information for non-DO or CC crosses. \emph{Not required for DO.}
	 
    When method = intensity:
       x: Numeric matrix, num.samples x num.snps, with X intensities for all founders and FALSE1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all founders and FALSE1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   code: Character vector containing two letter genotype codes for each founder sample. Sample IDs must be in names.
	   
    When method = allele:
       geno: Character matrix, num.samples x num.snps, with allele calls for all founders and FALSE1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   code: Character vector containing two letter genotype codes for each founder sample. Sample IDs must be in names.
  }
  \item{sex}{
    Character that is either FALSE or M, indicating the sex of the animals. Only used on X chromosome.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
  \item{chr}{
    Character indicating the current chromosome.
  }
  \item{trans.prob.fxn}{
    FALSEunction to use to determine the transition probability between markers.
  }
}
\author{
  Daniel Gatti
}
\keyword{ internal }
