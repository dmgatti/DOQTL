\name{create.Rdata.files}
\alias{create.Rdata.files}
\title{
  Convert *.txt files to *.Rdata files.
}
\description{
  This is used to convert the *.txt genotype probability files to *.Rdata files.
}
\usage{
  create.Rdata.files(prob.files)
}
\arguments{
  \item{prob.files}{
    Character vector containing genotype probability file names.
  }
}
\value{
  No value returned. Write out a *.Rdata file for each *.txt file.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ create.Rdata.files(prob.files) }
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
