\name{scanone.eqtl}
\alias{scanone.eqtl}
\alias{matrixeqtl.snps}
\title{
  Mapping using the Matrix EQTL algorithm.
}
\description{
  MatrixEQTL uses a series of matrix operations to greatly accelerate QTL mapping. It can accomodate additive covariates and a common kinship matrix for all phenotypes.
}
\usage{
  scanone.eqtl(expr, probs, K, addcovar, snps, sex)
}
\arguments{
  \item{expr}{
    Numeric matrix of phenotype values with samples in rows and phenotypes in columns. Rownames must contain sample IDs and colnames must contain phenotype names.
  }
  \item{probs}{
    Numeric three dimensional array containing the founder haplotype contributions. Num.samples by num.founders by num.markers.
  }
  \item{K}{
    Numeric matrix of kinship values for all samples. Num.samples by num.samples.
  }
  \item{addcovar}{
    Numeric matrix of additive covariates.
  }
  \item{snps}{
    Data.frame containing 4 columns with marker location information. SNP ID, chr, Mb, cM in columns 1 through 4, respectively.
  }
  \item{sex}{
    Character vector containing either "M" or F, indicating the sex of each sample. Used for mapping on the X chromosome.
  }
}
\details{
  Matrix EQTL centers and rotates the phenotype and genotype matrices using matrix operations. It only calculates the LOD score at each marker and does not provide coefficients.
}
\value{
  Numeric matrix of LOD scores for all phenotypes and markers.
}
\seealso{
  \code{\link{scanone}}
}
\references{
  Matrix eQTL: ultra fast eQTL analysis via large matrix operations.
  Shabalin AA.
  Bioinformatics. 2012 May 15;28(10):1353-8.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ scanone.eqtl(expr, probs, K, addcovar, snps) }
}
\keyword{ QTL }
