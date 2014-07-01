\name{qtl.qtlrel}
\alias{qtl.qtlrel}
\title{Use QTLRel to map a set of traits}
\description{
  This function accepts phenotypes, genotype probabilities and a sample kinship matrix and maps the requested traits using an eight state linear model.  FALSEixed covariates may be passed in as well.  The output is written to two files: *.LOD.txt (containing the LOD scores for each SNP) and *.coef.txt. (containing the model coefficient at each SNP).
}
\usage{
  qtl.qtlrel(pheno, probs, K, addcovar, intcovar, snps)
}
\arguments{
  \item{pheno}{
    Data frame, containing the sample IDs, phenotype data and covariates.
  }
  \item{probs}{
    3D numeric array, containing the genotype probabilities for all samples at each SNP.  Dimensions must be samples by states by SNPs and all dimensions must be named.
  }
  \item{K}{
    Numeric matrix, containing the kinship between individuals as computed by QTLRel.
  }
  \item{addcovar}{
    Numeric matrix containing additive covariates.
  }
  \item{intcovar}{
    Numeric matrix containing covariates that interact with the QTL.
  }
  \item{snps}{
     Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
}
\value{
  A list containing two elements:
  \item{lod}{
    Data.frame containing the SNP locations and LOD and p-values.
  }
  \item{coef}{
    Data.frame containing the model coefficients.
  }
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{plot.doqtl}}, \code{\link{scanone}}, \code{\link{scanone.perm}}
}
\examples{
  \dontrun{
    qtl.qtlrel(pheno, prob, K, covar = NULL, pheno.name = "")    
  }
}
\keyword{ QTL }
\keyword{ QTLRel }
\keyword{ mapping }
