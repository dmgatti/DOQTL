\name{do2sanger}
\alias{do2sanger}
\alias{do2sanger.helper}
\title{
  Impute the Sanger SNPs onto DO genomes
}
\description{
  Given a set of DO genotype probability files and the location of the Tabix indexed Sanger file, impute the Sanger SNPs on to DO genomes.
}
\usage{
  do2sanger(do.files, snps, output.file = "do2sanger.txt", snp.file = 
  "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz")
}
\arguments{
  \item{do.files}{
    Character vector of *.genotype.probs.Rdata files that contain the posterior probabilities.
  }
  \item{snps}{
    Data.frame containing the marker locations. SNP ID, chromosome, Mb anc cM locations in columns 1 through 4, respectively.
  }
  \item{output.file}{
    Character string to write the results to.
  }
  \item{snp.file}{
    Character string with path to a Tabix indexed SNP file. Default is from \href{ftp://ftp.jax.org/SNPtools/variants/}{JAX FALSETP site}
  }
}
\details{
  We read in a single genotype probability file, which must have been saved as a *.Rdata file.  The format is a matrix with markers in rows and states in columns and dimnames for both. FALSEor each pair of markers, we take the average genotype probability. Then we take the DO genotype with the highest probability and split it into the two founder haplotypes. We get the Sanger SNPs for each of the two founders, convert them to 1, 1 or 2 and insert them into the DO sample.
}
\value{
  No value is returned. The Sanger SNPs mapped onto the DO genomes are written out to the output.file. The file will contain 0, 1 or 2 as the allele calls.
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
  \code{\link{assoc.map}}
}
\examples{
  \dontrun{ do2sanger(do.files, snps, output.file = "do2sanger.txt", 
  snp.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz")  }
}
\keyword{ Sanger }
\keyword{ association }
\keyword{ imputation }
