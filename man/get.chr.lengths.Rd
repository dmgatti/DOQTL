\name{get.chr.lengths}
\alias{get.chr.lengths}
\title{
  Get chromosome lengths for the mouse
}
\description{
  Using the org.Mm.eg.db package, get the chromosome lengths.
}
\usage{
  get.chr.lengths()
}
\details{
  Keeps only the lengths of the chromosomes, not random or unmapped.
}
\value{
  Named vector of chromosome lengths.
}
\author{
  Daniel Gatti
}
\examples{
  get.chr.lengths()
}
\keyword{ chromosome }
