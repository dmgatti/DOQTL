\name{get.pattern.variants}
\alias{get.pattern.variants}
\title{get.pattern.variants}
\description{
  Given a vector of strains names and a set of varaints, return the variants for which the given strains have one allele and the remaining strains have the other allele.
}
\usage{
  get.pattern.variants(variants, strain.subset = NULL)
}
\arguments{
  \item{variants}{
    Matrix of data.frame of variants with four header columns containing annotation and snps in the remaining columns.
  }
  \item{strain.subset}{
    character vector of strain names that comprise one subset.
  }
}
\details{
  Use this function to obtain a subset of SNPs for which one set of strains contain one allele and the rest contain the other allele.  This might come up in the course of QTL mapping where several strains have a high allele and others have a low allele.
}
\value{
  Data.frame of variants that match the requested pattern.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{get.variants}}
}
\examples{
  data(example.snps)
  snp.subset = get.pattern.variants(example.snps[1:1000,], 
  strain.subset = c("A/J", "NOD/ShiLtJ", "NZO/HlLtJ", "AKR/J"))
}
