\name{add.sig.thr}
\alias{add.sig.thr}
\title{
  Add the significance thresholds to an existing QTL plot.
}
\description{
  Given a set of thresholds, add the autosomal and X chromosome thresholds to an existing QTL plot. 
  The user may call plot() and then this function, but \code{\link{plot.doqtl}} contains a 'sig.thr'
  argument and will call this automatically. The threhsolds can be obtained from
  \code{\link{get.sig.thr}}
}
\usage{
  add.sig.thr(sig.thr, sig.col = "red", chrsum)
}
\arguments{
  \item{sig.thr}{
    Numeric matrix or a numeric vector with the significance thredholds, typically obtained
    from \code{\link{get.sig.thr}}. If sig.thr is a matrix, then it must have 2 columns
    names 'A' and 'X' and each significance threshold is in one row.
  }
  \item{sig.col}{
    Numeric vector containing the plotting colors for each threshold. The length of sig.col
    must equal nrow(sig.thr).
  }
  \item{chrsum}{
    Numeric vector with the cumulative sum of the chromosome lengths. Must be named with
    chromosome names.
  }
}
\value{
  Plots autosomal and X chromosome thresholds on a QTL plot. Behavior on other plots is undetermined.
}
\seealso{
  \code{\link{get.sig.thr}}
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ 
     qtl = scanone(pheno = pheno, probs = probs, addcovar = addcovar, snps = anps)
     perms = scanone.perm(pheno = pheno, probs = probs, addcovar = addcovar, snps = anps)
     sig.thr = get.sig.thr(perms)
     plot(qtl)
     add.sig.thr(sig.thr, chrsum = cumsum(get.chr.lengths()))
     # Could also run:
     plot(qtl, sig.thr = sig.thr)
  }
}
\keyword{ permutations }
\keyword{ significance }
