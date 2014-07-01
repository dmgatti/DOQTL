\name{gene.plot}
\alias{gene.plot}
\alias{line.up.genes}
\alias{resolve.collisions}
\title{gene.plot}
\description{
  Given a genomic region, plot the genes in the interval. line.up.genes and resolve.collisions are internal functions.
}
\usage{
  gene.plot(mgi, col = "black", ...)
  line.up.genes(locs)
  resolve.collisions(locs) 
}
\arguments{
  \item{mgi}{
    Data.frame with MGI gene locations as returned by \code{\link{get.mgi.features}{get.mgi.features}}.
  }
  \item{col}{
    Color vector that is the color to use to plot genes. May be a single color for all genes or a vector with colors for each gene.
  }
  \item{locs}{
    Data.frame containing the location and names of genes. Columns are symbol, gene start Mb, gene end Mb, transcript start Mb, transcript end Mb and row.
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
