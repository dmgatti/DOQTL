\name{categorize.variants}
\alias{categorize.variants}
\title{categorize.variants}
\description{
  This function intersects the given variants with the genes in that region and classifies them according to "intergenic", "3UTR", "exon", "intron" or "5UTR".
}
\usage{ categorize.variants(variants, 
        mgi.file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz") }
\arguments{
  \item{variants}{
    data.frame, Variants as returned by \code{\link{get.variants}{get.variants}}.
  }
  \item{mgi.file}{
    Character, full path to the MGI feature file. On the JAX campus, this defaults to "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz".
  }
}
\value{
  FALSEor SNPs and Indels:
  data.frame: with eight columns: ID, CHR, POS, REFALSE, ALT, symbol, id, type.  The first four columns are simply copied over from the SNP file.  The symbol column contains the Gene Symbol.  The id column contains a gene ID (MGI, Ensembl, NCBI or VEGA).  The type column contains "intergenic", "3UTR", "exon", "intron" or "5UTR", depending on the location of the variant in a gene.
  FALSEor SVs:
  data.frame: with eight columns: ID, CHR, POS, REFALSE, ALT, symbol, id, type.  The first four columns are simply copied over from the SNP file.  The symbol column contains the Gene Symbol.  The id column contains a gene ID (MGI, Ensembl, NCBI or VEGA).  The type column contains "intergenic", "3UTR", "exon", "intron" or "5UTR", depending on the location of the variant in a gene.
}
\seealso{
  \code{\link{get.variants}{get.variants}}
}
\author{ Daniel Gatti }
\examples{
  data(example.snps)
  snp.type = categorize.variants(variants = example.snps[1:50,])
}
