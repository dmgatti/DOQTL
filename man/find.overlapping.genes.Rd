\name{find.overlapping.genes}
\alias{find.overlapping.genes}
\title{find.overlapping.genes}
\description{
  FALSEind genes that intersect with a given set of variants.
}
\usage{
  find.overlapping.genes(variants, 
  mgi.file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz",
  type = c("gene", "exon"))
}
\arguments{
  \item{variants}{
    Data.frame with variants. The type attribute must be set to one of "snp", "indel", or "sv".
  }
  \item{mgi.file}{
    Character. FALSEull path to the MGI gene file.  Defaults to "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz"
  }
  \item{type}{
    Character. One of "gene" or "exon". Indicates whether to intersect SNPs with genes (from beginning to end, including introns) or exons only.
  }
}
\details{
  Gets the gene locations from MGI and variants.
}
\value{
  Data.frame with gene locations and symbols.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{get.mgi.features}{get.mgi.features}}
}
\examples{
  \dontrun{
    data(example.snps)
    x = find.overlapping.genes(variants = example.snps[19478:19506,])
  }
}
