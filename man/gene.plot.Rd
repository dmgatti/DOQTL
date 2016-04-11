\name{gene.plot}
\alias{gene.plot}
\title{gene.plot}
\description{
  Given a genomic region, plot the genes in the interval. line.up.genes and resolve.collisions are internal functions.
}
\usage{
  gene.plot(mgi, rect.col = "grey30", text.col = "black", ...)
}
\arguments{
  \item{mgi}{
    Data.frame with MGI gene locations as returned by \code{\link{get.mgi.features}{get.mgi.features}}.
  }
  \item{rect.col}{
    Color vector that is the color to use to plot gene rectangles. May be a single color for all genes or a vector with colors for each gene. Default = "grey30".
  }
  \item{text.col}{
    Color vector that is the color to use to plot gene names. May be a single color for all genes or a vector with colors for each gene. Default = "black".
  }
  \item{...}{
    Other arguments to be passed into plot.
  }
}
\details{
  The spacing algorithm attempts to organize the genes in such a way that they do not collide. 
  The rect.col and text.col arguments are recycled if they are shorter than the number of genes. They can be used to highlight specific genes.
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
