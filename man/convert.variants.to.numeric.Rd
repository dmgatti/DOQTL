\name{convert.variants.to.numeric}
\alias{convert.variants.to.numeric}
\title{convert.variants.to.numeric}
\description{
Given a matrix or data frame with character {ACGT} SNP values, convert them to numeric values with the major allele coded as 0 and the minor allele as 1.
}
\usage{
convert.variants.to.numeric(variants)
}
\arguments{
  \item{variants}{
     Data frame with four header columns. The SNPs must be in {ACGT} format in columns 5 through ncol(snps).
  }
}
\details{
This function is used before calling plot.alleles() to convert the SNPs into a numeric form suitable for plotting.
}
\value{
Data.frame of the same dimentions as the snps argument, but with the alleles converted to 0 or 1.
}
\author{
Daniel Gatti
}
\examples{
  data(example.snps)
  numeric.snps = convert.variants.to.numeric(example.snps[1:100,])
}
