\name{batch.normalize}
\alias{batch.normalize}
\alias{normalize.batches}
\alias{quantilenorm}
\title{Batch normalize the X & Y intensity data.}
\description{
  This function batch normalizes the X & Y intensity data by subtracting batch medians from the X & Y intensities.
}
\usage{
  batch.normalize(path = ".", snps)
  quantilenorm(x1, y1, x2, y2)
}
\arguments{
  \item{path}{
     Character, the full path to the input files, which must be either "x.txt" and "y.txt" or "x.filt.txt" and "y.filt.txt".
  }
  \item{snps}{
     Data.frame, with three columns containing SNP ID, chromosome and Mb location in that order. May be obtained from \href{ftp://ftp.jax.org/MUGA}{ftp://ftp.jax.org/MUGA}.
  }
  \item{x1}{
    Numeric matrix containing X intensities for batch 1 containing samples in rows and markers in columns. Number of samples should be larger than x2.
  }
  \item{y1}{
    Numeric matrix containing Y intensities for batch 1 containing samples in rows and markers in columns. Number of samples should be larger than y2.
  }
  \item{x2}{
    Numeric matrix containing X intensities for batch 2 containing samples in rows and markers in columns.
  }
  \item{y2}{
    Numeric matrix containing Y intensities for batch 2 containing samples in rows and markers in columns.
  }
}
\details{
  quantile.norm adjusts the intensities of samples in batch 2 to those of batch 1. The number of samples in batch 1 should be greater than the number of samples in batch 2. At each SNP, we form quantiles of the X1 (or Y1) intensity distribution, discarding the upper and lower 0.01% to reduce the influence of outliers. We then adjust the X2 (or Y2) intensities to match the X1 (or Y1) quantiles. This will not work well if batch 2 has too few samples. Empirically, it seems that a minimum of 50 samples are needed.
}
\value{
  FALSEor batch.normalize: returns value is returned.  The batch normalized intensities are written to "x.filt.batch.norm.txt" and "y.filt.batch.norm.txt".
  
  FALSEor quantilenorm: returns normalized X and Y values for batch 2.
}
\author{
  Daniel Gatti
}
\note{
  FALSEuture releases may include more sophisticated normalization algorithms.
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
