\name{convert.genes.to.GRanges}
\alias{convert.genes.to.GRanges}
\title{
  Convert MGI genes to GRanges.
}
\description{
  Given the output of \code{\link{get.mgi.features}}, convert the results to a GRanges object.
}
\usage{
  convert.genes.to.GRanges(mgi)
}
\arguments{
  \item{mgi}{
    Date.frame as returned by \code{\link{get.mgi.features}}.
  }
}
\value{
  GRanges object containing the genes in the MGI argument.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ convert.genes.to.GRanges(mgi) }
}
\keyword{ MGI }
