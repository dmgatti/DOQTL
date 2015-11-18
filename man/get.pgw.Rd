\name{get.pgw}
\alias{get.pgw}
\title{
  Get the genome wide p-value.
}
\description{
  Given a set of permutations and a LOD or -log10(p-value), return the adjusted genom-wide p-value. This requires autosomal and X chromosome permutation values as provided by \link{\code{scanone.perm}}.
If you are using p-values, make sure to convert both the permutations and the p-values to -log10(-p-value).
}
\usage{
  get.pgw(stat, chr, perms)
}
\arguments{
  \item{stat}{
    Numeric vector of mapping statistics. Must be either LOD or -log10(p-values).
  }
  \item{chr}{
    Numeric vector of same legnth as stat containing the chromosome on which each statistic occurs.
  }
  \item{perms}{
    Numeric matrix with the maximum LOD (or maximum -log10(p-value)) from each permuation.
  }
}
\value{
  Numeric vector containing the adjusted genome-wide p-value for each statistic.
}
\author{
  Daniel Gatti, Petr Simecek
}
\seealso{
  \code{\link{get.sig.thr}}
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
