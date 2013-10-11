\name{snp.plot}
\alias{snp.plot}
\title{snp.plot}
\description{
  Lower level function that makes a tile plot the given SNPs with the major allele colored dark and the minor allele light.  Optionally, add SNPs that match a certain pattern and genes and a QTL score.
}
\usage{
  snp.plot(variants, col = c("black", "grey50", "white"), cluster = T, ref, highlight, 
  pattern.snps, mgi, qtl)
}
\arguments{
  \item{variants}{
    Data.frame with variants as returned by \code{\link{get.variants}}.
  }
  \item{col}{
    Color vector with SNP colors for no call, alternate allele and reference allele.
  }
  \item{cluster}{
    Boolean that indicates if the strains should be clustered.
  }
  \item{ref}{
    Chracter that is the reference strain to use.  Must be present in the strains in variants.
  }
  \item{highlight}{
    Character vector with strain names to highlight in the plot. Strain names must be present in the strains in variants.
  }
  \item{pattern.snps}{
    Data.frame with SNPs that match some pattern. These will be plotted below the SNPs.
  }
  \item{mgi}{
    Data.frame with genes in the interval to be plotted.
  }
  \item{qtl}{
    data.frame with QTL values containing the chromosome, bp position and QTL score in column 1 to 3. Default = NULL.
  }
}
\details{
  Different strains can be used as the reference strain by using the ref argument. Otherwise, the major allele is plotted dark and the minor allele lighter.  The QTL values will be scaled within the plotting interval and drawn with black (low) and red (high). If a strain pattern is provided, the SNPs matching the pattern are plotted in orange and any genes that they intersect with are also colored orange. Otherwise, genes are colored blue.
}
\value{
  Produces a tile plot with the locations along the horizontal axis and the strains in the SNP matrix along the vertical axis.  Also, if a set of strains is given in the pattern argument, the SNPs that match that pattern are returned, categorized according to which gene they lie within.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{variant.plot}}, \code{\link{get.mgi.features}}, \code{\link{categorize.variants}}
}
\examples{
  \dontrun{
  data(qtl)
  strains = get.strains()
  variants = get.variants(chr = 7, start = 103, end = 105, 
  strains = strains[c(2,4,8,10,15:18)])
  pattern.snps = get.pattern.variants(variants, strain.subset = c("C57BL/6J", "NOD/ShiLtJ", 
  "NZO/HlLtJ"))
  mgi = get.mgi.features(chr = 7, start = 103, end = 105, source = "MGI", type = "gene")
  variants = convert.variants.to.numeric(variants)
  snp.plot(variants = variants, pattern.snps = pattern.snps, mgi = mgi, qtl = qtl)
  }
}
