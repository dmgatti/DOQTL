\name{gene.plot}
\alias{gene.plot}
\title{gene.plot}
\description{
  Given a genomic region, plot the genes in the interval.
}
\usage{
  gene.plot(mgi, col = "black", ...)
}
\arguments{
  \item{mgi}{
    Data.frame with MGI gene locations as returned by \code{\link{get.mgi.features}{get.mgi.features}}.
  }
  \item{col}{
    Color vector that is the color to use to plot genes. May be a single color for all genes or a vector with colors for each gene.
  }
  \item{...}{
    Other arguments to be passed into plot.
  }
}
\details{
  The spacing algorithm attempts to organize the genes in such a way that they do not collide. Currently, this works 90% of the time. Resizing the plot resolves this sometimes.
}
\value{
  Data.frame with gene locations and symbols.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{get.mgi.features}{get.mgi.features}}
}
\examples{
  \dontrun{
  data(example.genes)
  g = gene.plot(mgi = example.genes)
  }
}
