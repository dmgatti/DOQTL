\name{filter.samples}
\alias{filter.samples}
\title{Filter X, Y and genotype data by call rate}
\description{
  This function reads in “x.txt”, “y.txt”, “geno.txt” and “call.rate.batch.txt” from the user specified input directory and removes samples with a call rate less than the supplied threshold (default = 0.9).  It then writes the files out to “x.filt.txt”, “y.filt.txt”, “geno.filt.txt” and “call.rate.batch.filt.txt”.
}
\usage{
  filter.samples(path = ".", thr = 0.9)
}
\arguments{
  \item{path}{
    Character, the full path to the directory where the files reside.
  }
  \item{thr}{
    Numeric, call rate threshold below which samples will be removed. Default = 0.9.
  }
}

\value{
  data.frame with the sample IDs and call rates of the removed samples.  The intensity and genotype files are written to the input directory.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{extract.raw.data}}, \code{\link{batch.normalize}}
}
\examples{
  filter.samples(path = "/tmpdir/output")
}
\keyword{ MUGA }
\keyword{ call rate }
\keyword{ normalization }
