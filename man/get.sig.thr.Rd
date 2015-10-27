\name{get.sig.thr}
\alias{get.sig.thr}
\title{
  Get the significance thresholds.
}
\description{
  Given a set of permutations, return the significance thresholds on the autosomes and X chromosome.
}
\usage{
  get.sig.thr(perms, alpha = 0.05, Xchr = TRUE)
}
\arguments{
  \item{perms}{
    Numeric matrix with the maximum LOD (or maximum -log10(p-value)) from each permuation.
  }
  \item{alpha}{
    Numeric vector containing the alpha level of the significance thresholds.
  }
  \item{Xchr}{
    Logical value (default = TRUE) that indicates that the permutations contain
    boht autosomal and X chromosome thresholds.
  }
}
\value{
  Numeric matrix containing the significance thresholds for the autosomes
  in column 1 an dthe X chromossome in column 2.
}
\author{
  Daniel Gatti
}
\references{
  The X chromosome in quantitative trait locus mapping.
  Broman KW, Sen S, Owens SE, Manichaikul A, Southard-Smith EM, Churchill GA.
  Genetics. 2006 Dec;174(4):2151-8.
}
\examples{
  \dontrun{ 
    perms = scanone.perm(pheno = pheno, probs = probs, addcovar = addcovar, snps = anps)
    get.sig.thr(perms) 
  }
}
\keyword{ permutations }
\keyword{ significance }
