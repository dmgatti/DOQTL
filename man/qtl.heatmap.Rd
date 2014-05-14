\name{qtl.heatmap}
\alias{qtl.heatmap}
\title{Plot a Heatmap of all QTL}
\description{
  This function accepts LOD scores and SNP locations and plots a heatmap of the QTL.  Each QTL is scaled to between 0 and 1 and the QTL are hierarchically clustered.  A tile plot is then created with QTL colored from white (strong) to black (weak).  The motivation for this plot is to visualize QTL that may be coincident across different phenotypes.
}
\usage{
  qtl.heatmap(lod, ...)
}
\arguments{
  \item{lod}{
    Data frame, containing the SNP names, chromosome an dbp locations in columns 1 through 3 and LOD scores for all phenotyped in the remaining columns.
  }
  \item{...}{
    Additional arguments to be passed to plot.
  }
}
\value{
  No value is returned.  A tile plot is plotted from the data.
}
\author{
  Daniel Gatti
}
\note{
  The QTL the appear in white are *not* neccessarily signficant.  They are simply the strongest QTL for each phenotype.
}
\seealso{
  \code{\link{plot.doqtl}}, \code{\link{scanone}}, \code{\link{scanone.perm}}
}
\examples{
  \dontrun{
    qtl.heatmap(lod)
  }
}
\keyword{ QTL }
\keyword{ QTLRel }
\keyword{ mapping }
