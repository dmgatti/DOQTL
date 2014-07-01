\name{get.gene.name}
\alias{get.gene.name}
\title{
  Get the gene symbol
}
\description{
  Given an MGI ID, get the gene symbol from an MGI data.frame.
}
\usage{
  get.gene.name(value, mgi)
}
\arguments{
  \item{value}{
    Character string containing an MGI ID.
  }
  \item{mgi}{
    Data.frame containing MGI data from \code{\link{get.mgi.features}}.
  }
}
\value{
  Character string containing the gene symbol.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{get.mgi.features}}
}
\examples{
  \dontrun{ get.gene.name(value, mgi) }
}
\keyword{ MGI }
