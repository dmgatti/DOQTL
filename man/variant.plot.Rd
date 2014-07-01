\name{variant.plot}
\alias{variant.plot}
\title{variant.plot}
\description{
  Convenience function that makes a tile plot the given SNPs with the major allele colored dark and the minor allele light.  Optionally, add SNPs that match a certain pattern and genes and a QTL score. This calls \code{\link{sdp.plot}} after gathering all of the data.
}
\usage{
  variant.plot(var.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz", 
  mgi.file = "ftp://ftp.jax.org/SNPtools/genes/MGI.20130703.sorted.txt.gz", chr, 
  start, end, type = c("snp", "indel", "sv"), strains = c("A/J", "C57BL/6J", 
  "129S1/SvImJ", "CAST/EiJ", "NOD/ShiLtJ", "NZO/HlLtJ", "PWK/PhJ", "WSB/EiJ"), 
  ref = "C57BL/6J", pattern, qtl)
}
\arguments{
  \item{var.file}{
    Character string with the full path to the tabix indexed Sanger/UNC SNP file. Default is a file at the Center for Genome Dynmaics at the Jackson Laboratory.
  }
  \item{mgi.file}{
    Character string with the full path to the tabix indexed MGI feature file. Default is a file at the Center for Genome Dynmaics at the Jackson Laboratory.
  }
  \item{chr}{
    Character with the chromosome ID to plot on.  Required.
  }
  \item{start}{
    Numeric value with the start position to plot at. May be in bp or Mb. If the value is 200 or less, it is assumed to be in Mb.  Otherwise it is assumed to be in bp. Required.
  }
  \item{end}{
    Numeric value with the end position to plot at. May be in bp or Mb. If the value is 200 or less, it is assumed to be in Mb.  Otherwise it is assumed to be in bp. Required.
  }
  \item{type}{
    Character indicating the type of variant to retrieve. One of "snp", "indel", or "sv". Default = "snp".
  }
  \item{strains}{
    Character vector with the strains to use.  Must be from the set of strains returned by \code{\link{get.strains}{get.strains}}. Default = Collaborative Cross/Diversity Outbred founders.
  }
  \item{ref}{
    Character string with the reference strain. Must be one of the names in the strains vector argument. Default = C57BL/6J.
  }
  \item{pattern}{
    Character vector with strain names in the strains argument. We search for SNPs for which these strains have one allele and all other strains have the opposite allele. A more sophisticated pattern matching method may be included in future releases. Default = NULL.
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
  \code{\link{convert.variants.to.numeric}}, \code{\link{get.mgi.features}}, \code{\link{categorize.variants}}, \code{\link{snp.plot}}
}
\examples{
  \dontrun{
  data(qtl)
  strains = get.strains()[c(2,4,8,10,15:18)]
  variant.plot(chr = 7, start = 103, end = 105, strains = strains, 
  pattern = c("C57BL/6J", "NOD/ShiLtJ", "NZO/HlLtJ"), qtl = qtl)
  }
}
