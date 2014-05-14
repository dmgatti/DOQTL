\name{intensity.plots}
\alias{founder.F1.intensity.plot}
\alias{intensity.mean.covar.plot}
\title{
  Plot founders and F1 hybrids or genotype state means and variances on an intensity plot.
}
\description{
  Given the X and Y (or rho and theta) array intensities, plot the values for one SNP and color and mark the founders and FALSE1s.
}
\usage{
  founder.F1.intensity.plot(theta, rho, s = 1, is.founder.F1, plotNew = TRUE, ...)
  intensity.mean.covar.plot(s, states, theta, rho, r.t.means, r.t.covars, sample, ...)
}
\arguments{
  \item{theta}{
    Numeric matrix containing theta (or X) intensity values.
  }
  \item{rho}{
    Numeric matrix containing rho (or Y) intensity values.
  }
  \item{s}{
    Numeric integer containing the SNP index in rho and theta to plot.
  }
  \item{is.founder.F1}{
    Boolean vector containing TRUE for samples that are DO founders or FALSE1s and FALSEALSE for those that are not.
  }
  \item{plotNew}{
    Boolean that is TRUE if a new plot should be made and FALSEALSE if the points should be added to the existing plot.
  }
  \item{states}{
    Character vector containing genotypes state letter codes.
  }
  \item{r.t.means}{
    Numeric three dimensional array containing genotype state means in rho and theta coordinates.
  }
  \item{r.t.covars}{
    Numeric three dimensional array containing genotype state variances in rho and theta coordinates.
  }
  \item{sample}{
    Numeric value containing the sample index into the rho and theta matrices to plot in orange.
  }
  \item{\dots}{
    Additional arguments to be passed to plot.
  }
}
\value{
  No return value. A scatter plot is made of the rho and theta (or X and Y) intensities at the requested SNP. Founders and F1 hybrid samples are colored.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{ founder.F1.intensity.plot(theta, rho, s = 1, is.founder.F1, plotNew = TRUE, ...) }
}
\keyword{ intensities }
