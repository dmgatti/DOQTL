\name{pxg.plot}
\alias{pxg.plot}
\title{
  Phenotype by genotype plot at a single marker.
}
\description{
  This function plots the phenotype versus the most probable genotype at a single marker.
}
\usage{
  pxg.plot(pheno, pheno.col, probs, snp.id, snps, legend = TRUE, sex = NA, covar)
}
\arguments{
  \item{pheno}{
    Data.frame containing the phenotype data with samples in rows and phenotypes in columns. Sample IDs in rownames and phenotype names in colnames.
  }
  \item{pheno.col}{
    Numeric or character vector: Either a vector of number that indicate columns to use or a set of column names in pheno.
  }
  \item{probs}{
    3D numeric array, containing the founder haplotype contributions or genotype probabilities. The sample IDs, founder letter codes and markers IDs must be in dimnames.
  }
  \item{snp.id}{
    Character string containing the marker to plot at.  Must be in SNPs and dimnames(probs)[[3]].
  }
  \item{snps}{
    Data.frame containing 4 columns with marker location information. SNP ID, chr, Mb, cM in columns 1 through 4, respectively.
  }
  \item{legend}{
    Boolean that is true if a legend explaining the founder letter codes should be plotted.
  }
  \item{sex}{
    Character that is either FALSE or M, indicating sex. This is used only when the SNP ID is on the X chromosome and all of the samples are male.  In this case, there are only 8 genotype states.
  }
  \item{covar}{
    Vector of categories for each sample (i.e. diet, treatement, etc.) which will be plotted as different symbols. Optional.
  }
}
\details{
  The most probable genotype is inferred from the eight founder allelic contributions at the marker. If a founder has a value > 0.75, it is assumed to be homozygous.
}
\value{
  Creates a plot with the phenotype on the y-axis and the 36 DO genotypes on the x-axis.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{plot.doqtl}}, \code{\link{coefplot}}
}
\examples{
  \dontrun{ pxg.plot(pheno, pheno.col, founder.probs, snp.id, snps) }
}
\keyword{ QTL }
