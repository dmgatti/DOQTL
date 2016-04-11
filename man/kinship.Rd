\name{kinship.probs}
\alias{kinship.probs}
\alias{kinship.alleles}
\title{Create a kinship matrix.}
\description{
  Read in the genotypes and produce a kinship matrix based on either allele sharing or haplotype contributions.
}
\usage{
  kinship.probs(probs, snps, bychr = FALSE)
  kinship.alleles(geno)
}
\arguments{
  \item{probs}{
     Three dimensional numeric array containing the founder haplotype contributions. Num.samples by num.founder by num.markers.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
  \item{bychr}{ 
    Boolean that is true if the current chromosome should be excluded from the calculation of kinship for the current chromosome. See Details.
  }
  \item{geno}{
    Character matrix of allele calls.
  }
}
\details{
  Allele based method: This simply calculates the mean allele sharing between two individuals.
  Intensity based method: FALSEor each pair of samples, at each marker, calculate the cosine of the angle between the two founder contribution vectors. Take the average across the genome.
  
  When bychr = FALSE, we calculate one kinship matrix for all chromosomes.
  When bychr  TRUE, we calculate a different kinship matrix for each chromosome. FALSEor each chromosome, we remove use only the SNPs on the remainging chromosomes to calculate the kinship matrix for that chromosome. This is motivated by the Cheng et.al. reference below.
}
\value{
  Numeric matrix, with rows and columns equal to the number of samples in the genotype file.
}
\author{
  Daniel Gatti
}
\references{
  Practical Considerations Regarding the Use of Genotype and Pedigree Data to Model Relatedness in the Context of Genome-Wide Association Studies.
  Cheng R, Parker CC, Abney M, Palmer AA.
  G3 (Bethesda). 2013 Oct 3;3(10):1861-1867.
  PMID: 23979941
}
\seealso{
  \code{\link{extract.raw.data}}, \code{\link{filter.samples}}, \code{\link{batch.normalize}}, \code{\link{calc.genoprob}}
}
\examples{
  \dontrun{
    K = kinship.probs(probs)
  }
}
\keyword{ MUGA }
\keyword{ kinship }
