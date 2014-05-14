\name{example.pheno}
\alias{example.pheno}
\docType{data}
\title{
  Example phenotypes.
}
\description{
  Example phenotypes.
}
\usage{ example.pheno }
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
  head(example.pheno)
}
\keyword{datasets}
