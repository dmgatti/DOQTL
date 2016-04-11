\name{scanone}
\alias{scanone}
\title{Perform a genome scan.}
\description{
  This function is the main QTL mapping function for point mapping at each marker. The user must supply phenotypes, genotype probabilities and marker locations. Optional kinship, additive and interactive covariates may be passed in. Scanone regresses the phenotype on the genotype probabilities and produces a LOD score and founder allele effects at each marker.
}
\usage{
  scanone(pheno, pheno.col = 1, probs, K, addcovar, intcovar, snps, model = 
  c("additive", "full"))
}
\arguments{
  \item{pheno}{
    data.frame containing phenotype data. Required. rownames must contain sample IDs and there must be a column labelled 'sex' to perform correct mapping on the X chromosome.
  }
  \item{pheno.col}{
    numeric or character vector: Either a vector of number that indicate columns to use or a set of column names in pheno.
  }
  \item{probs}{
    Numeric three dimensional array, containing the founder haplotype contributions or genotype probabilities. The sample IDs, founder letter codes and markers IDs must be in dimnames.
  }
  \item{K}{
    numeric matrix, containing the additive kinship matrix. The samples IDs must be in rownames and colnames.
  }
  \item{addcovar}{
    data.frame or numeric matrix, containing any additive covariates. Sample IDs must be in rownames.
  }
  \item{intcovar}{
    data.frame or numeric matrix, containing any covariates that interact with the QTL. Sample IDs must be in rownames.
  }
  \item{snps}{
    data.frame containing 4 columns with marker location information. SNP ID, chr, Mb, cM in columns 1 through 4, respectively.
  }
  \item{model}{
    character string, containing one of "additive" or "full", indicating the type of model to fit.
  }
}
\value{
  A list containing two elements:
  \item{lod}{
     Data.frame with nine columns containing the marker locations, LOD scores and p-values for each phenotype. Column names are SNP_ID, Chr, Mb_NCBI38, cM, perc.var, lrs, lod, p and neg.log10.p.
  }
  \item{coef}{
     Numeric matrix containing the founder allele effects at each locus. Colnames are the additive covariates plus the founder terms.
  }
}
\details{
  We require a sex argument in addcovar for mapping on the X chromosome. This is true even if all of your samples are of the same sex. It also means that mapping on the autosomes will occur with sex as an additive covariate. See A Guide to QTL Mapping with R/qtl by Karl Broman, pages 108-118.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{scanone.perm}}
}
\examples{
  \dontrun{scanone(pheno, pheno.col = 1, probs)}
}
\keyword{ QTL }
