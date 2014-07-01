\name{get.variants}
\alias{get.variants}
\title{get.variants}
\description{
  Get SNPs for a set of strains for a specific regions of chromosomes.
}
\usage{
  get.variants(file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz", 
  chr, start, end, type = c("snp", "indel", "sv"), strains, polymorphic  =  TRUE, 
  quality)
}
\arguments{
  \item{file}{
    Character, full path to the SNP file to use.  Default is the file at the Center for Genome Dynamics at The Jackson Laboratory.
  }
  \item{chr}{
    Numeric vector, chr for each start and end position.  Chr, start and end must all have the same length.
  }
  \item{start}{
    Numeric vector, start position in Mb or bp for each chr.  Chr, start and end must all have the same length.
  }
  \item{end}{
    Numeric vector, end position in Mb or bp for each chr.  Must be greater than or equal to the corresponding start value.  Chr, start and end must all have the same length.
  }
  \item{type}{
    Character, with one of "snp", "indel", "sv". Indicates the type of variantes being queried.
  }
  \item{strains}{
    Character vector, listing the strains to retrieve. The names can be obtained from get.strains(). Default returns all strains.
  }
  \item{polymorphic}{
    Boolean. Default  =  TRUE. If TRUE, retrieve only SNPs that are polymorphic among the requested strains.  If FALSEALSE, return all SNPs in the requested interval(s). 
  }
  \item{quality}{
    Integer, denoting the confidence levels to return. This looks at the quality column of all requested strains and takes values greater than or equal to the given number. Look at the documentation for the SNP file you are using for details on the quality score. Default returns all SNPs.
  }
}
\value{
  A list of data.frames for each region requested.  Each data.frame will contain five header columns; ID, CHROM, POS, REFALSE & ALT corresponding to the SNP ID, chromosome, bp position, reference and alternate alleles.  There follows two columns for each requested strain containing the allele calls and quality scores for each strain. If there is only one requested region, a single data frame is returned.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{get.strains}}
}
\examples{
  \dontrun{
  available.strains = get.strains()
  strains = available.strains[c(2, 5, 8, 9, 13:15, 17)]
  snps = get.variants(chr = 7, start = 103, end = 105, strains = strains)
  }
}
