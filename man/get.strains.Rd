\name{get.strains}
\alias{get.strains}
\title{get.strains}
\description{
  Get available strain names for a set of strains from a variant file.
}
\usage{
  get.strains(
  file = "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz")
}
\arguments{
  \item{file}{
    Character, full path to the variants file to use.  Default is the file at the Center for Genome Dynamics at The Jackson Laboratory.
  }
}
\value{
  Character vector of strain names in the variant file.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{
  available.strains = get.strains()
  strains = available.strains[c(2, 5, 8, 9, 13:15, 17)]
  snps = get.variants(chr = 7, start = 103, end = 105, strains = strains)
  }
}
