\name{estimate.cluster.params}
\alias{estimate.cluster.params}
\alias{keep.homozygotes}
\title{
  Estimate genotype cluster means and variances
}
\description{
  Given the X and Y intensity data, perform model based clustering and estimate the genotype state cluster means and variances.
}
\usage{
  estimate.cluster.params(founders, data, chr)
  keep.homozygotes(founders)
}
\arguments{
  \item{founders}{
    List containing founder information for non-DO or CC crosses. \emph{Not required for DO.}
	 
    When method = intensity:
       x: Numeric matrix, num.samples x num.snps, with X intensities for all founders and FALSE1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all founders and FALSE1s (if available). Sample IDs and SNP IDs must be in rownames and colnames.
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
    Character containing the current chromosome.
  }
}
\details{
  At each marker, use mclust to perform model based clustering on all of the data and get estimates of the means and variances for each cluster. Then assign each of the 36 genotype states to the nearest founder cluster. 
  
  keep.homozygotes is an internal helper function.
}
\value{
   List containing two elements:
   \item{r.t.means}{
	 Three dimensional array containing rho and theta genotype cluster means.
   }
   \item{r.t.covars}{
	 Three dimensional array containing rho and theta genotype cluster variances.
   }
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{hmm.intensity}}
}
\examples{
  \dontrun{ estimate.cluster.params(founders, data, chr) }
}
\keyword{ cluster }
