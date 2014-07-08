\name{scanone.perm}
\alias{scanone.perm}
\title{Perform a genome scan.}
\description{
  This function is the main QTL mapping function for point mapping at each marker. The user must supply phenotypes, genotype probabilities and marker locations. Optional kinship, additive and interactive covariates may be passed in. scanone.perm regresses the phenotype on the genotype probabilities and produces a LOD score and founder allele effects at each marker.
}
\usage{
  scanone.perm(pheno, pheno.col = 1, probs, addcovar, intcovar, snps, model = 
  c("additive", "full"), path = ".", nperm = 1000)
}
\arguments{
  \item{pheno}{
    data.frame containing phenotype data. Required. rownames must contain sample IDs and there must be a column labelled 'sex' to perform correct mapping on the X chromosome.
  }
  \item{pheno.col}{
    numeric or character vector: Either a vector of number that indicate columns to use or a set of column names in pheno.
  }
  \item{probs}{
    3D numeric array, containing the founder haplotype contributions or genotype probabilities. The sample IDs, founder letter codes and markers IDs must be in dimnames.
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
    Character string, containing one of "additive" or "full", indicating the type of model to fit.
  }
  \item{path}{
    Character string containing the location where permutation results should be written.
  }
  \item{nperm}{
    Numeric, containing the number of permutations to run.
  }
}
\value{
  A list containing two elements:
  \item{lod}{
     Data.frame with seven columns containing the marker locations and LOD scores for each phenotype.
  }
  \item{coef}{
     Numeric matrix containing the founder allele effects at each locus.
  }
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{scanone.perm}}
}
\examples{
  \dontrun{scanone.perm(pheno, pheno.col = 1, probs)}
}
\keyword{ QTL }
