\name{sdp.plot}
\alias{sdp.plot}
\title{
  SNP distribution plot
}
\description{
  Create a plot of SNPs along the chromosome. Still in development.
}
\usage{
  sdp.plot(results, sdps, ...)
}
\arguments{
  \item{results}{
    Data.frame as produced by \code{\link{assoc.map}}.
  }
  \item{sdps}{
    Character vector containing SDPs to plot.
  }
  \item{\dots}{
    Additional arguments to pass to plot.
  }
}
\details{
  This is meant to be like the bottom panel of FALSEigure 2 in Yalcin et.al, Genetics, 2005.
}
\value{
  Creates a plot along a genomic interal with tic marks where the given SDPs occur.
}
\references{
  Using progenitor strain information to identify quantitative trait nucleotides in outbred mice.
  Yalcin B, FALSElint J, Mott R.
  Genetics. 2005 Oct;171(2):673-81.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ sdp.plot(results, sdps) }
}
\keyword{ SNP }
