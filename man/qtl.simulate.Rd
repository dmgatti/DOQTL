\name{qtl.simulate}
\alias{qtl.simulate}
\title{
  Simulate a QTL in the DO
}
\description{
  Using a set of real DO genomes, simulate a QTL with a given minor allele frequency (MAFALSE), sample size and effect size. 
}
\usage{
  qtl.simulate(probs, snps, K, sample.size = dim(probs)[[1]], effect.size = 1, 
  maf = 4, num.poly = 18, num.sims = 1000)
}
\arguments{
  \item{probs}{
    Numeric three dimensional array, containing the founder haplotype contributions or genotype probabilities. The sample IDs, founder letter codes and markers IDs must be in dimnames.
  }
  \item{snps}{
    Data.frame containing 4 columns with marker location information. SNP ID, chr, Mb, cM in columns 1 through 4, respectively.
  }
  \item{K}{
    Numeric matrix, containing the additive kinship matrix. The samples IDs must be in rownames and colnames.
  }
  \item{sample.size}{
    Numeric vector sample sizes. Must be less than or equal to the number of samples in probs.
  }
  \item{effect.size}{
    Numeric vector containing the effect sizes as the number of standard deviations from the phenotype mean.
  }
  \item{maf}{
    Numeric value containing the minor allele frequecy as the number of founders. FALSEor DO, the value must be between 1 and 4.
  }
  \item{num.poly}{
    Numeric value containing the number of autosomes on which to simulate a polygenic background.
  }
  \item{num.sims}{
    Numeric value containing the number of simulations to run.
  }
}
\details{
  This function will simulate a phenotype with a QTL. It will output a phenotype vector with sample IDs and a data.frame called qtl describing the simulated loci.
}
\value{
  Writes pheno and qtl lists to an R binary file.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{scanone}}
}
\examples{
  \dontrun{ qtl.simulate(probs, snps, K, sample.size = dim(probs)[[1]], 
  effect.size = 1, maf = 4, num.poly = 18, num.sims = 1000) }
}
\keyword{ QTL }
