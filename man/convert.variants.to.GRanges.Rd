\name{convert.variants.to.GRanges}
\alias{convert.variants.to.GRanges}
\title{convert.variants.to.GRanges}
\description{
  Given a data.frame of SNPs, convert the SNP locations to a GRanges object.
}
\usage{
  convert.variants.to.GRanges(variants)
}

\arguments{
  \item{variants}{
     Data frame with four header columns and SNPs in the remaining columns.
  }
}
\details{
  This function creates a GRanges object from the CHR and POS columns of the SNP data.frame.
}
\value{
  GRanges object containing the SNP locations.
}
\author{
Daniel Gatti
}
\examples{
  data(example.snps)
  gr = convert.variants.to.GRanges(example.snps[1:100,])
}

