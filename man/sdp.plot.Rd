\name{sdp.plot}
\alias{sdp.plot}
\title{
  Plot association mapping results.
}
\description{
  Plot the founder strain distribution patterns (SDP) of the SNPs in results.
}
\usage{
  sdp.plot(results, ...)
}
\arguments{
  \item{results}{
    Data.frame containing output from \code{\link{assoc.map}}.
  }
  \item{\dots}{
    Additional arguments passed to plot.
  }
}
\details{
  Given the output from \code{\link{assoc.map}}, plot the SDPs for the SNPs in results. This will show which strains have the minor allele at each SNP in results. You may also add this plot to the top of \code{\link{assoc.plot}}.
}
\value{
  Produces a plot. There is no return value.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{assoc.map}}, \code{\link{assoc.plot}}
}
\examples{
  \dontrun{ 
    results = assoc.map(pheno = pheno, pheno.col = 1, probs = probs, K = K, addcovar = addcovar, 
	snps = snps, chr = 1, start = 40, end = 45)
    sdp.plot(results[results[,12] > 3,])
  }
}
\keyword{ association }
\keyword{ qtl }
\keyword{ sdp }
