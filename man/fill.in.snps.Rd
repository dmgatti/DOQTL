\name{fill.in.snps}
\alias{fill.in.snps}
\title{
  Interpolate between SNPs at the same cM value.
}
\description{
  Go through the SNPs and look for stretches where the differences in cM values from one SNP to the next equals 0. Interpolate from the first SNP to the last SNP in such a stretch.
}
\usage{
  fill.in.snps(snps)
}
\arguments{
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
}
\details{
  On each chromosome, look for sets of contiguous markers that have the same cM location. Interpolate evenly spaced markers spanning the markers proximal and distal to this set.
}
\value{
  Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ fill.in.snps(snps) }
}
\keyword{ snps }
