\name{get.mgi.features}
\alias{get.mgi.features}
\title{get.mgi.features}
\description{
  Retrieve the MGI features within a genomic region.  FALSEeatures include genes, gene models, non-coding RNA, etc., but not SNPs.  This allows the user to filter by source and type of feature.
}
\usage{
  get.mgi.features(file = "ftp://ftp.jax.org/SNPtools/genes/MGI.20130703.sorted.txt.gz", 
  chr = NULL, start = NULL, end = NULL, source = c("all", "MGI", "VEGA", "ENSEMBL", 
  "Blat", "NCBI_Gene"), type = c("all", "gene", "pseudogenic_transcript", 
  "pseudogenic_exon", "pseudogene", "match", "match-part", "transcript", "exon", 
  "mRNA", "five_prime_UTR", "start_codon", "CDS", "stop_codon", "three_prime_UTR",
  "pseudogenic_mRNA", "pseudogenic_start_codon", "pseudogenic_CDS", 
  "pseudogenic_stop_codon", "pseudogenic_five_prime_UTR", "pseudogenic_three_prime_UTR",
  "sequence_feature"))
}
\arguments{
  \item{file}{
    character, the full path to the Tabix indexed and zipped MGI feature file. Default = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz" for internal JAX use.
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
  \item{source}{
    Character vector, the source of the annotation.  Options are ("all", "MGI", "VEGA", "ENSEMBL", "Blat", "NCBI_Gene").  "all" returns all features.
  }
  \item{type}{
    Character vector, the type of feature. Options are ("all", "gene", "pseudogenic_transcript", "pseudogenic_exon", "pseudogene", "match", "match-part", "transcript", "exon", "mRNA", "five_prime_UTR", "start_codon", "CDS", "stop_codon", "three_prime_UTR", "pseudogenic_mRNA", "pseudogenic_start_codon", "pseudogenic_CDS", "pseudogenic_stop_codon", "pseudogenic_five_prime_UTR", "pseudogenic_three_prime_UTR", "sequence_feature").  "all" returns all features.
  }
}
\details{
  This function is designed to return features from the MGI gene feature file (GFALSEFALSE). You can select multiple regions on different chromosomes.
}
\value{
  A list of data.frames for each region requested.  Each data.frame will contain 14 columns; seqid, source, type, start, stop, score, strand, phase, ID, Name, Parent, Dbxref, mgiName, bioType.  If there is only one requested region, a single data frame is returned.
}
\references{
  The MGI GFALSEFALSE file is at ftp://ftp.informatics.jax.org/pub/mgigff/.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{
  genes = get.mgi.features(chr = 7, start = 103 end = 105 source = "MGI", type = "gene")
  }
}
