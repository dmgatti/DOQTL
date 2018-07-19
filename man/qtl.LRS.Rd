\name{qtl.LRS}
\alias{qtl.LRS}
\alias{permutations.qtl.LRS}
\title{
  QTL mapping with no kinship.
}
\description{
  These are two functions that perform QTL mapping without a kinship matrix. They use the QR decomposition to speed up the computation.  Other than for a quick screen or for assessing significance thresholds, we do not recommend mapping without a kinship matrix. They are included for historical reasons
}
\usage{
  qtl.LRS(pheno, probs, snps, addcovar = NULL)
  permutations.qtl.LRS(pheno, probs, snps, addcovar, nperm = 1000, 
  return.val = c("lod", "p"))
}
\arguments{
  \item{pheno}{
    Data.frame containing the phenotype data with samples in rows and phenotypes in columns. Sample IDs in rownames and phenotype names in colnames.
  }
  \item{probs}{
    Numeric three dimenaional array, containing the founder haplotype contributions or genotype probabilities. The sample IDs, founder letter codes and markers IDs must be in dimnames.
  }
  \item{snps}{
    data.frame containing 4 columns with marker location information. SNP ID, chr, Mb, cM in columns 1 through 4, respectively.
  }
  \item{addcovar}{
    data.frame or numeric matrix, containing any additive covariates. Sample IDs must be in rownames.
  }
  \item{nperm}{
    Numeric value containing the number of permutations to run.
  }
  \item{return.val}{
    Character string containing either "LRS" or "p", indicating the type of return statistic.
  }
}
\details{
  The function performs Haley-Knott regression at the markers using the founder haplotype contributions in probs.
}
\value{
  List containing two elements. LRS: a data.frame containing the LOD scores. coef: numeric matrix containing model coefficients.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{scanone}}
}
\examples{
  \dontrun{ qtl.LRS(pheno, probs, snps, addcovar = NULL) }
}
\keyword{ QTL }
