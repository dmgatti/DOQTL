\name{fast.qtlrel}
\alias{fast.qtlrel}
\title{
  QTL mapping using QTLRel
}
\description{
  This extracts teh core of the QTLRel algorithm for additive covariates.
}
\usage{
  fast.qtlrel(pheno, probs, K, addcovar, snps)
}
\arguments{
  \item{pheno}{
    Data.frame containing the phenotype data. Sample IDs in rownames.
  }
  \item{probs}{
    Three dimensional numeric array containing the founder halplotype contributions. Num.samples by num.founder by num.markers.  Sample IDs, founder letters and SNP IDs must be in dimnames.
  }
  \item{K}{
    Numeric matrix containing the kinship between samlpes. Sample IDs must be in rownames and colnames.
  }
  \item{addcovar}{
    Numeric matrix containing additive covariates to run in the mapping model. Sample IDs must be in rownames.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
}
\details{
  We extracted code from QTLRel, but removed several options to speed up the pipeline for QTL mapping with additive covariates and a kinship matrix.
}
\value{
  List containing two elements:
  \item{lod}{
     List containing two elements:
	 A: Numeric matrix containing the SNP ID, chromosome, Mb, cM, percent variance explained, likelihood ratio statistic, log of the odd ratio, p-value and -log10(p-value) for autosomes.
	 X: Numeric matrix containing the SNP ID, chromosome, Mb, cM, percent variance explained, likelihood ratio statistic, log of the odd ratio, p-value and -log10(p-value) for X chromosome.
  }
  \item{coef}{
    List containing two elements:
	A : Numeric matrix containing QTL model coefficients on autosomes. Markers in rows.
	X : List containing two elements for the X chromosome:
	  FALSE : Numeric matrix containing QTL model coefficients for females. Markers in rows.
	  M : Numeric matrix containing QTL model coefficients for males. Markers in rows.
  }
}
\references{
  Cheng R, Abney M, Palmer AA, Skol AD.  
  QTLRel: an R package for genome-wide association studies in which relatedness is a concern.
  BMC Genet. 2011 Jul 27;12:66. 
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{scanone}}, \code{\link{scanone}}
}
\examples{
  \dontrun{ fast.qtlrel(pheno, probs, K, addcovar, snps) }
}
\keyword{ QTL }
