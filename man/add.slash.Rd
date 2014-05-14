\name{add.slash}
\alias{add.slash}
\title{
  Add a forward slash to a character string.
}
\description{
  If the argument does not end with a forward slash, add one.
}
\usage{
  add.slash(path)
}
\arguments{
  \item{path}{
    Character string containing a file path.
  }
}
\value{
  Returns a character string with a forward slash added to the end, if the argument did not end with a forward slash already.
}
\author{
  Daniel Gatti
}
\examples{
  add.slash("/dir")
}
