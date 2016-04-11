\name{add.missing.F1s}
\alias{add.missing.F1s}
\title{Add Missing F1 Samples}
\description{
  Given a set of CC or DO founders, impute the genotypes or intensities of missing F1s.
}
\usage{
  add.missing.F1s(founders, snps, sampletype = c("DO", "CC", "DOF1", "HS", "other"))
}
\arguments{
  \item{founders}{
    List, required: Contains either an element called 'geno' or two elements called 'x' and 'y'.
  }
  \item{snps}{
    Data.frame, required: Data.frame containing the SNPs. SNP ID, chr, Mb and cM locations in columna 1 - 4, respectively.
  }
  \item{sampletype}{
    Character string, required: indicates the type of crss. One of "DO", "CC", "DOF1", "HS", "other".
  }
}
\value{
  List with the founders data structure updated to include missing F1 samples.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{
    load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
    founders = add.missing.F1s(founders, snps = muga_snps)
  }
}
