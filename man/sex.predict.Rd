\name{sex.predict}
\alias{sex.predict}
\title{Determine the sex of each sample}
\description{
  This function uses the mean X and Y chromosome intensities to predict the sex of each sample.  Linear discriminant analysis from the MASS package is used.
}
\usage{
  sex.predict(x, y, snps, plot = FALSE)
}
\arguments{
  \item{x}{
     Numeric matrix, num.samples x num.snps, containing X intensities for all samples.
  }
  \item{y}{
     Numeric matrix, num.samples x num.snps, containing Y intensities for all samples.
  }
  \item{snps}{
     Data.frame containing SNP IDs, chromosomes, Mb and cM locations in columns 1 through 4, respectively.
  }
  \item{plot}{
     Boolean that will create a plot of mean X and Y Chr intensities if TRUE.  Default = FALSEALSE.
  }
}
\value{
  Character vector with sex assignments.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{
    data(snps)
    sex.predict(x, y, sex = rep(F, ncol(x)))
  }
}
\keyword{ sex }
