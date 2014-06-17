\name{plot.doqtl}
\alias{plot.doqtl}
\title{Plot a QTL}
\description{
  Given a genome scan produced by scanone, create a plot of the LOD score across all chromosomes.
}
\usage{
  \method{plot}{doqtl}(x, stat.name = c("lod", "neg.log10.p"),  sig.thr = NULL, sig.col = "red", ...)
}
\arguments{
  \item{x}{
    DOQTL object containing the LOD score and model coefficients.
  }
  \item{stat.name}{
    Character string containing the name of the mapping statistic.
  }
  \item{sig.thr}{
    Numeric vector containing significance thresholds.
  }
  \item{sig.col}{
    Colors to use when plotting the significance thresholds above. 
  }
  \item{...}{
    Additional arguments to pass to plot.
  }
}
\value{
  Creates a QTL plot.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{
    qtl = scanone(pheno = pheno, pheno.col = 1, probs = probs, snps = snps)
	plot(qtl)
  }
}
