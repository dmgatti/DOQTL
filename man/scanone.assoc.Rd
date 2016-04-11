\name{scanone.assoc}
\alias{scanone.assoc}
\title{
  Map the entire genome using association mapping. 
}
\description{
  Use the imputed Sanger SNPs and the DO haplotype probabilities. Between each 
  pair of markers, get the unique founder strain distribution patters (SDPs)
  and use the DO haplotype probabilities to compute the DO SNPs.
}
\usage{
  scanone.assoc(pheno, pheno.col, probs, K, addcovar, markers, cross = c("DO", "CC", "HS"), sdp.file, ncl)
}
\arguments{
  \item{pheno}{
    Data.frame containing the phenotype values in columns and samples in rows. rownames(pheno) must contain the sample IDs.
  }
  \item{pheno.col}{
    numeric or character vector: Either a vector of number that indicate columns to use or a set of column names in pheno.
  }
  \item{probs}{
    Numeric three dimensional array containing the founder haplotype contributions. Num.samples by num.founders by num.markers.
  }
  \item{K}{
    List containing numeric matrices of leave-one-chromosome-out kinship values for all samples. Num.samples by num.samples.
  }
  \item{addcovar}{
    Numeric matrix of additive covariates.
  }
  \item{cross}{
    Character string indicating the type of cross. One of "DO", "CC", or "HS".
  }
  \item{markers}{
    Data.frame containing 4 columns with marker location information. SNP ID, chr, Mb, cM in columns 1 through 4, respectively.
  }
  \item{sdp.file}{
    Character string containing the full path to the condensed SDP file. This file is created using \code{\link{condense.sanger.snps}}.
  }
  \item{ncl}{
    Integer containing the number of cores to use for parallel execution.
  }
}
\details{
  This function imputes the Sanger SNPs onto DO genomes and performs association 
  mapping. The speed relies upon a support file that is created using
  \code{\link{condense.sanger.snps}}. The support file is a tab-delimited,
  Tabix indexed file that contains the chromosome, bp position and SDP for
  each polymorphic Sanger SNP. It treats tri-morphic SNPs and heterozygotes
  as alternate alleles.
}
\value{
  GRanges object containing the p-value for each SNP.
}
\seealso{
  \code{\link{scanone}}
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ scanone.assoc(pheno, pheno.col, probs, K, addcovar, markers, cross = c("DO", "CC", "HS"), sdp.file, ncl) }
}
\keyword{ QTL }
\keyword{ scanone }
