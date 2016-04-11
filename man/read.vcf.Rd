\name{read.vcf}
\alias{read.vcf}
\alias{get.vcf.strains}
\title{
  Read and parse VCF data
}
\description{
  *** CURRENTLY UNDER CONSTRUCTION ***. This function is designed to read in the Sanger VCF files for SNPs, Indels and structural variants. 
}
\usage{
  read.vcf(vcf.file, gr, chr = 1, start = 4, end = 4.5, strains,
           return.val = c("allele", "number"), return.qual = TRUE, csq = FALSE)
}
\arguments{
  \item{vcf.file}{
    Character string containing full path to a Sanger SNP, Indel or structural variant VCF file.
  }
  \item{gr}{
    GRanges object containing one or more ranges to query. Use either this argument or the chr, start and end arguments, but not both.
  }
  \item{chr}{
    Character or numeric vector containing mouse chromosome IDs to query.
  }
  \item{start}{
    Numeric vector containing the start position on each chromosome to query. Must be the same length as chr. If the value is less than 200, it is assumed to be in Mb. Values greater than 200 are assumed to be in bp.
  }
  \item{end}{
    Numeric vector containing the end position on each chromosome to query. Must be the same length as chr. If the value is less than 200, it is assumed to be in Mb. Values greater than 200 are assumed to be in bp.
  }
  \item{strains}{
    Character vector containing strain names, as listed in the vcf.file, to retrieve. If missing, then all strains are returned. Use \code{\link{get.vcf.strains}} to get the strain names from the VCF file.
  }
  \item{return.val}{
    Character string that is either "allele", which returns the nucleotide allele calls, or "number", which returns numeric values.
  }
  \item{return.qual}{
    Boolean that is TRUE if the quality columns should be returned.
  }
  \item{csq}{
    Boolean that is TRUE if the consequence column should be returned.
  }
}
\details{
  There is a very nice package called vcf2geno, but it is designed to work with the 1000 Genomes VCF format. The Sanger format is slightly different, hence the need for this function. Future plans include the addition of variant consequence selection.
  The Sanger variant files have 'snp', 'indel' or 'SV' in them and we use this to determine the type of file being processed.
}
\value{
  If one chromosome location is requested, a data.frame containing the requested variants and quality scores. The variant locations and reference alleles are in the first 5 columns and the variants and quality scores for each strain are in the remaining columns. If more than one location is requested, a list of data.frames containing the variants.
}
\references{
  The \href{http://www.sanger.ac.uk/resources/mouse/genomes/}{Sanger Mouse Genomes Project} generated these variants.
  Raw data is at the \href{ftp://ftp-mouse.sanger.ac.uk/}{Sanger FALSETP site}.
 
    Next-generation sequencing of experimental mouse strains.
    Yalcin B, Adams DJ, FALSElint J and Keane TM
    Mammalian genome : official journal of the International Mammalian Genome Society 2012;23;9-10;490-8
    PUBMED: 22772437
	
    Sequencing and characterization of the FALSEVB/NJ mouse genome.
    Wong K, Bumpstead S, Van Der Weyden L, Reinholdt LG, Wilming LG, Adams DJ and Keane TM
    Genome biology 2012;13;8;R72
    PUBMED: 22916792
	
    The fine-scale architecture of structural variants in 17 mouse genomes.
    Yalcin B, Wong K, Bhomra A, Goodson M, Keane TM, Adams DJ and FALSElint J
    Genome biology 2012;13;3;R18
    PUBMED: 22439878
	
    The genomic landscape shaped by selection on transposable elements across 18 mouse strains.
    Nellaker C, Keane TM, Yalcin B, Wong K, Agam A, Belgard TG, FALSElint J, Adams DJ, FALSErankel WN and Ponting CP
    Genome biology 2012;13;6;R45
    PUBMED: 22703977
	
    High levels of RNA-editing site conservation amongst 15 laboratory mouse strains.
    Danecek P, Nellaker C, McIntyre RE, Buendia-Buendia JE, Bumpstead S, Ponting CP, FALSElint J, Durbin R, Keane TM and Adams DJ
    Genome biology 2012;13;4;26
    PUBMED: 22524474
	
    Mouse genomic variation and its effect on phenotypes and gene regulation.
    Keane TM, Goodstadt L, Danecek P, White MA, Wong K, Yalcin B, Heger A, Agam A, Slater G, Goodson M, FALSEurlotte NA, Eskin E, Nellaker C, Whitley H, Cleak J, Janowitz D, Hernandez-Pliego P, Edwards A, Belgard TG, Oliver PL, McIntyre RE, Bhomra A, Nicod J, Gan X, Yuan W, van der Weyden L, Steward CA, Bala S, Stalker J, Mott R, Durbin R, Jackson IJ, Czechanski A, Guerra-Assuncao JA, Donahue LR, Reinholdt LG, Payseur BA, Ponting CP, Birney E, FALSElint J and Adams DJ
    Nature 2011;477;7364;289-94
    PUBMED: 21921910
	
    Sequence-based characterization of structural variation in the mouse genome.
    Yalcin B, Wong K, Agam A, Goodson M, Keane TM, Gan X, Nellaker C, Goodstadt L, Nicod J, Bhomra A, Hernandez-Pliego P, Whitley H, Cleak J, Dutton R, Janowitz D, Mott R, Adams DJ and FALSElint J
    Nature 2011;477;7364;326-9
    PUBMED: 21921916
 
}
\author{
  Daniel Gatti
}
\note{
  !!! This function is still under development !!!
}
\examples{
  \dontrun{ 
    read.vcf(vcf.file = 
    "ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v3.snps.rsIDdbSNPv137.vcf.gz",
    chr = 1, start = 5, end = 5.5)
    read.vcf(vcf.file = 
    "ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v3.snps.rsIDdbSNPv137.vcf.gz",
    gr = GRanges(seqnames = 1, ranges = IRanges(start = 5, end = 5.5)))
  }
}
\keyword{ Sanger }
\keyword{ VCF }
\keyword{ variants }
