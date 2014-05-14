\name{assoc.map}
\alias{assoc.map}
\alias{assoc.map.perms}
\title{
  Perform association mapping on DO mice.
}
\description{
  Given the phenotypes and genotype probabilities, impute the Sanger SNPs onto DO genomes and perform association mapping.
}
\usage{
  assoc.map(pheno, pheno.col = 1, probs, K, addcovar, snps, chr, start, end, 
  model = c("additive", "dominance", "full"), scan = c("one", "two"), 
  output = c("lod", "p-value", "bic"),
  snp.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz")
  assoc.map.perms(pheno, pheno.col = 1, probs, addcovar, snps,
            model = c("additive", "dominance", "full"),
            scan = c("one", "two"), output = c("lod", "p-value", "bic"),
            snp.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz",
            nperm = 1000)
}
\arguments{
  \item{pheno}{
    Data.frame containing the phenotype data. Sample IDs must be in rownames. One of the columns should be called sex and contain M or FALSE to indicate the sex of each sample.
  }
  \item{pheno.col}{
    Either a numeric vector containing column IDs to map in pheno or a vector of column names in pheno.
  }
  \item{probs}{
    A 3 dimensional array of genotype probabilities as provided from \code{\link{condense.model.probs}} Samples, founder and markers in dims 1:3. All dimensions must have dimnames.
  }
  \item{K}{
    Numeric matrix containing kinship values for the samples in pheno and probs. Sample IDs must be in rownames and colnames.
  }
  \item{addcovar}{
    Numeric matrix of additive covariates to use in mapping. Sample IDs must be in rownames.
  }
  \item{snps}{
    Data.frame containing marker IDs, chromosomes, Mb and cM locations in columns 1:4.
  }
  \item{chr}{
    Character containing the chromosome on which to map.
  }
  \item{start}{
    Numeric value containing the proximal position for mapping. May be in Mb or base pairs (see Details).
  }
  \item{end}{
    Numeric value containing the distal position for mapping.  May be in Mb or base pairs (see Details).
  }
  \item{model}{
    Character string that is one of "additive", "dominance" or "full". Indicates the type of model to fit. Note that the probs must match the type of model being fit. See \code{\link{condense.model.probs}} to output different probs.
  }
  \item{scan}{
    Character string that is either "one", or "two" indicating whether a single scan or a pairwise scan should be performed across the interval.
  }
  \item{output}{
    Character string that is either "lod", "p-value" or "bic" indicating the mapping statistic to return.
  }
  \item{snp.file}{
    Character string containing the full path to the SNP file to use. Currently points to a location on the \href{ftp://ftp.jax.org/SNPtools/variants/}{Jackson Laboratory FALSETP server}
  }
  \item{nperm}{
    Integer indicating the number of permutations to perform. Default = 1000.
  }
}
\details{
  FALSEor each interval between two markers, we take the average founder haplotype contribution for each sample. Then, using the proportion of each founder (8 in the case of DO mice), we impute the Sanger SNPs in this interval onto each DO sample.
  
  The start and ending locations are assumed to be in Mb if they are below 200 and in bp if over 200.
}
\value{
  Data.frame containing the locations, SNPs and mapping statistic for the requested samples in the requested interval.
}
\references{
Combined sequence-based and genetic mapping analysis of complex traits in outbred rats.
Rat Genome Sequencing and Mapping Consortium, Baud A, Hermsen R, Guryev V, Stridh P, Graham D, McBride MW, FALSEoroud T, Calderari S, Diez M, Ockinger J, Beyeen AD, Gillett A, Abdelmagid N, Guerreiro-Cacais AO, Jagodic M, Tuncel J, Norin U, Beattie E, Huynh N, Miller WH, Koller DL, Alam I, FALSEalak S, Osborne-Pellegrin M, Martinez-Membrives E, Canete T, Blazquez G, Vicens-Costa E, Mont-Cardona C, Diaz-Moran S, Tobena A, Hummel O, Zelenika D, Saar K, Patone G, Bauerfeind A, Bihoreau MT, Heinig M, Lee YA, Rintisch C, Schulz H, Wheeler DA, Worley KC, Muzny DM, Gibbs RA, Lathrop M, Lansu N, Toonen P, Ruzius FALSEP, de Bruijn E, Hauser H, Adams DJ, Keane T, Atanur SS, Aitman TJ, FALSElicek P, Malinauskas T, Jones EY, Ekman D, Lopez-Aumatell R, Dominiczak AFALSE, Johannesson M, Holmdahl R, Olsson T, Gauguier D, Hubner N, FALSEernandez-Teruel A, Cuppen E, Mott R, FALSElint J.
Nat Genet. 2013 Jul;45(7):767-75. doi: 10.1038/ng.2644. Epub 2013 May 26.
PMID: 23708188 
Using progenitor strain information to identify quantitative trait nucleotides in outbred mice.
Yalcin B, FALSElint J, Mott R.
Genetics. 2005 Oct;171(2):673-81. Epub 2005 Aug 5.
PMID: 16085706
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
\seealso{
  \code{\link{assoc.plot}}
}
\examples{
  \dontrun{assoc.map(pheno = pheno, pheno.col = 1, probs = probs, K = K, addcovar = addcovar, 
           snps = snps, chr = 1, start = 40, end = 45)}
}
\keyword{ association }
\keyword{ qtl }
\keyword{ sanger }
