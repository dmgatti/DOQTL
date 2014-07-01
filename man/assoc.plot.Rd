\name{assoc.plot}
\alias{assoc.plot}
\title{
  Plot association mapping results.
}
\description{
  After performing association mapping using \code{\link{assoc.map}}, plot the mapping statistic along with genes in the QTL interval.
}
\usage{
  assoc.plot(results, 
  mgi.file = "ftp://ftp.jax.org/SNPtools/genes/MGI.20130703.sorted.txt.gz", 
  highlight, highlight.col = "red", thr, ...)
}
\arguments{
  \item{results}{
    Data.frame containing output from \code{\link{assoc.map}}.
  }
  \item{mgi.file}{
    Character string containing the full path to a Tabix indexed file of gene locations. Default points to a version of the MGI genome feature file.
  }
  \item{highlight}{
    Character vector containing gene symbols to highlight in the plot.
  }
  \item{highlight.col}{
    Vector of colors to use when highlighting genes.
  }
  \item{thr}{
    Numeric value above which data points should be colored red and SNPs with these points returned.
  }
  \item{\dots}{
    Additional arguments passed to plot.
  }
}
\details{
  Given the output from \code{\link{assoc.map}}, plot the LOD or difference in BIC values across the QTL interval in the top panel. Plot the genes in the interval in the lower panel. 
}
\value{
  A plot with the mapping statistic in the top panel and genes in the lower panel.
  If thr is not missing, then filter the SNPs in the results argument and return only those with a mapping statistic greater than thr.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{assoc.map}}
}
\examples{
  \dontrun{ 
    results = assoc.map(pheno = pheno, pheno.col = 1, probs = probs, K = K, addcovar = addcovar, 
	snps = snps, chr = 1, start = 40, end = 45)
    assoc.plot(results)
  }
}
\keyword{ association }
\keyword{ qtl }
