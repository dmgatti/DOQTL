\name{call.genotypes}
\alias{call.genotypes}
\title{Call the most probably genotype by finding the Viterbi path through the genotype states for each sample.}
\description{
  This function is run after the HMM has computed genotype probabilities and produced genotype mean and covariances for each SNP.  It takes the HMM parameters and calculates the Viterbi path through the hidden genotype states.  The output is the most probable genotype call for each sample at each SNP.
}
\usage{
  call.genotypes(cross, theta.means, rho.means, theta.rho.covars)
}
\arguments{
  \item{x}{
       x: Numeric matrix, num.samples x num.snps, with X intensities for all samples.
  }
  \item{y}{
       y: Numeric matrix, num.samples x num.snps, with Y intensities for all samples.
  }
  \item{sex.gen}{
       sex.gen: data.frame, with two columns names 'sex' and 'gen'. The sex column must contain "M" or "F" and the gen column must contain numbers indicating the DO outbreeding generation for each sample.  Rownames must contain the sample IDs.
  }
  \item{theta.means}{
    numeric matrix, containing the theta genotypes state means. num.states rows by num.snps columns.
  }
  \item{rho.means}{
    numeric matrix, containing the rho genotypes state means. num.states rows by num.snps columns.
  }
  \item{theta.rho.covars}{
    numeric 3D array, containing the covariance matrices for each genotype cluster at each SNP.
  }
}

\value{
  character matrix, containing the Viterbi path through the genotype states for each sample at all SNPs. num.samples by num.snps.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{calc.genoprob}}
}
\examples{
  \dontrun{
    call.genotypes(cross, theta.means, rho.means, theta.rho.covars)
  }
}
\keyword{ MUGA }
\keyword{ HMM }
\keyword{ genotypes }
