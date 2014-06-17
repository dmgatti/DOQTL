\name{emission.probs.allele}
\alias{emission.probs.allele}
\alias{emission.probs.intensity}
\title{
  Calculate the emission probabilities 
}
\description{
  Calculate the emission probabilities 
}
\usage{
  emission.probs.allele(founders, chr, snps, sex = c("F", "M"))
  emission.probs.intensity(data, params)
}
\arguments{
  \item{founders}{
     List containing founder information for non-DO or CC crosses.
    When method = allele:
       geno: Character matrix, num.samples x num.snps, with allele calls for all founders and FALSE1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   code: Character vector containing two letter genotype codes for each founder sample. Sample IDs must be in names.
  }
  \item{data}{
    A list with named elements containing the information needed to reconstruct genomes. 
	
    When method = intensity:
       x: Numeric matrix, num.samples x num.snps, with X intensities for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all samples. Sample IDs and SNP IDs must be in rownames and colnames.
       sex: Character vector, containing "M" or F indicating sex. Sample IDs must be in names.
	   gen: Character matrix containing the generation of DO outbreeding for each sample. FALSEor the DO, this should be "DO" followed by a number with no space between them. FALSEor CC mice, this should be CC. Sample IDs must be in names.
  }
  \item{chr}{
    Character string containing the current chromosome.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
  \item{sex}{
    Character containing either FALSE or M, indicating the sex of the samples.
  }
  \item{params}{
    List containing two elements:
	r.t.means : Three dimensional array containing rho and theta genotype cluster means.
    r.t.covars : Three dimensional array containing rho and theta genotype cluster variances.
  }
}
\details{
  Given the founder allele calls or the rho and theta intensities, calculate the emission probabilities.
}
\value{
  Three dimensional numeric array containing the emission probabilities for each sample at each marker.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ emission.probs.allele(founders, chr, snps, sex = c(F, "M")) }
}
\keyword{ internal }
