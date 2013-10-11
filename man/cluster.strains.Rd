\name{cluster.strains}
\alias{cluster.strains}
\title{cluster.strains}
\description{
  Given a set of numeric SNPs, cluster the strains based on allele sharing. This function hierarchically clusters the strains based on the proportion of SNPs that share alleles between each strain.  Note: numeric snps, not character snps, are the input to this function.
}
\usage{cluster.strains(variants) }

\arguments{
  \item{variants}{
    data.frame, numeric variants as returned by \code{\link{convert.variants.to.numeric}}.
  }
}
\value{
  data.frame: of numeric variants with the strains clustered.
}
\seealso{
\code{\link{convert.variants.to.numeric}}, \code{\link{variant.plot}}
}
\author{ Daniel Gatti }
\examples{
  data(example.snps)
  variants = convert.variants.to.numeric(variants = example.snps[1:100,])
  variants = cluster.strains(variants)
}
