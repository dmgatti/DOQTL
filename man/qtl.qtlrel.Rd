\name{qtl.qtlrel}
\alias{qtl.qtlrel}
\title{Use QTLRel to map a set of traits}
\description{
  This function accepts phenotypes, genotype probabilities and a sample kinship matrix and maps the requested traits using an eight state linear model.  Fixed covariates may be passed in as well.  The output is written to two files: *.LOD.txt (containing the LOD scores for each SNP) and *.coef.txt. (containing the model coefficient at each SNP).
}
\usage{
  qtl.qtlrel(pheno, prob, K, covar = NULL, pheno.name = "", test = "None")
}
\arguments{
  \item{pheno}{
    Data frame, containing the sample IDs, phenotype data and covariates.
  }
  \item{prob}{
    3D numeric array, containing the genotype probabilities for all samples at each SNP.  Dimensions must be samples by states by SNPs and all dimensions must be named.
  }
  \item{K}{
    Numeric matrix, containing the kinship between individuals as computed by QTLRel.
  }
  \item{pheno.name}{
    Character, the phenotype name. This is appended as a prefix to the output files.
  }
  \item{test}{
    Character, default = "None", which performs a likelihood ratio test.  Other options include "F" or "Chisq".  See the QTLRel scanOne documentation for more details.
  }
}
\value{
  No value is returned.  The LOD scores and model coeffients are written to *.LOD.txt and *.coef.txt.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{qtl.plot}}, \code{\link{scanone}}, \code{\link{scanone.perm}}
}
\examples{
  \dontrun{
    qtl.qtlrel(pheno, prob, K, covar = NULL, pheno.name = "", test = "None")    
  }
}
\keyword{ QTL }
\keyword{ QTLRel }
\keyword{ mapping }

