\name{batch.normalize}
\alias{batch.normalize}
\title{Batch normalize the X & Y intensity data.}
\description{
  This function batch normalizes the X & Y intensity data by subtracting batch medians from the X & Y intensities.
}
\usage{
  batch.normalize(path = ".", snps)
}
\arguments{
  \item{path}{
     Character, the full path to the input files, which must be either "x.txt" and "y.txt" or "x.filt.txt" and "y.filt.txt".
  }
  \item{snps}{
     Data.frame, with three columns containing SNP ID, chromosome and Mb location in that order. May be obtained from \href{ftp://ftp.jax.org/MUGA}{ftp://ftp.jax.org/MUGA}.
  }
}
\value{
  No value is returned.  The batch normalized intensities are written to “x.filt.batch.norm.txt” and “y.filt.batch.norm.txt”.
}
\author{
  Daniel Gatti
}
\note{
  Future releases may include more sophisticated normalization algorithms.
}
\seealso{
  \code{\link{extract.raw.data}}, \code{\link{filter.samples}}
}
\examples{
  \dontrun{
    load(url("ftp://ftp.jax.org.MUGA/muga_snps.Rdata"))
    batch.normalize(path = "/demo/MUGA/", snps = muga_snps)
  }
}
\keyword{ MUGA }
\keyword{ intensity }
\keyword{ normalization }
