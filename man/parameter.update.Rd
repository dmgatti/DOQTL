\name{parameter.update}
\alias{parameter.update.alleles}
\alias{parameter.update.intensity}
\title{
  Parameter updating in HMM
}
\description{
  Not intended for external use.
}
\usage{
  parameter.update.alleles(geno, b, pseudocounts, prsmth)
  parameter.update.intensity(data, params, prsmth, founder.means)
}
\arguments{
  \item{geno}{
    Data.frame containing allele calls.
  }
  \item{b}{
    Three dimensional numeric array containing emission probabilities.
  }
  \item{pseudocounts}{
    Three dimensional numeric array containing pseudocounts for updating.
  }
  \item{prsmth}{
    Three dimensional numeric array containing posterior genotype probabilities.
  }
  \item{data}{
    A list with named elements containing the information needed to reconstruct genomes. 
	
    When method = intensity:
       x: Numeric matrix, num.samples x num.snps, with X intensities for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   gen: Character matrix containing the generation of DO outbreeding for each sample. FALSEor the DO, this should be "DO" followed by a number with no space between them. FALSEor CC mice, this should be CC. Sample IDs must be in names.
  }
  \item{params}{
    List containing two elements:
	r.t.means: three dimensional numeric array containing the genotype cluster means at each marker.
	r.t.covars: three dimensional numeric array containing the genotype cluster variances at each marker.
  }
  \item{founder.means}{
    Numeric matrix containing the founder intensity means for each marker.
  }
}
\author{
  Daniel Gatti
}
\keyword{ internal }
