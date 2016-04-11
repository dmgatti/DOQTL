\name{example.qtl}
\alias{example.qtl}
\docType{data}
\title{
  Example QTL.
}
\description{
  Example QTL for the albino and HDLD2 traits in \code{\link{example.pheno}}.  Albino is a binary, Mendelian trait that maps to the Tyrosinase locus on Chr 7. HDLD2 (high density lipoprotein) is a complex trait with many loci.
}
\usage{ example.qtl }
\format{
  A data frame with 149 observations on the following 8 variables. Sample IDs in rownames and phenotype names in colnames.
  \describe{
    \item{\code{Sample}}{
	  Character vector containing sample IDs.
	}
    \item{\code{Sex}}{
	  Character vector containing the sex of each animal.
	}
    \item{\code{Gen}}{
	   Character vector containing DO outbreeding generation and litter.
	}
    \item{\code{Diet}}{
	  Character vector containing diet, either chow or hf for high fat.
	}
    \item{\code{Coat.Color}}{
	  Character vector containing text description of mouse colors.
	}
    \item{\code{albino}}{
	  Numeric vector contining 1 if the mouse's coat color was white.
	}
    \item{\code{black}}{
	  Numeric vector contining 1 if the mouse's coat color was black.
	}
	\item{\code{HDLD2}}{
	  Numeric vector of high density lipoprotein values. 
	}
  }
}
\details{
  Data from Svenson et.al. paper below.
}
\references{
  High-resolution genetic mapping using the Mouse Diversity outbred population.
  Svenson KL, Gatti DM, Valdar W, Welsh CE, Cheng R, Chesler EJ, Palmer AA, McMillan L, Churchill GA.
  Genetics. 2012 FALSEeb;190(2):437-47
  PMID: 223445611
}
\examples{
  names(example.qtl)
  names(example.qtl[[1]])
}
\keyword{datasets}
