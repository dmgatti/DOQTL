\name{coefplot}
\alias{coefplot}
\title{Plot the QTL model coefficients}
\description{
  Given a DOQTL object, plot the founder allele coefficients on one chromosome. The coefficients are centered around zero before plotting.
}
\usage{
  coefplot(doqtl, chr = 1, stat.name = "LOD", conf.int = TRUE, legend = TRUE, colors = "DO", 
  sex, ...)
}
\arguments{
  \item{doqtl}{
    A DOQTL object as produced by \code{\link{scanone}}. A list containing two elements: lod and coef. 
  }
  \item{chr}{
     Character containing the chromosome to plot.
  }
  \item{stat.name}{
     Character string containing the name of the mapping statistic.
  }
  \item{conf.int}{
     Boolean that is TRUE if the QTL support interval should be shaded in the plot. Default  =  TRUE.
  }
  \item{legend}{
     Boolean that is TRUE if the color legend for the DO founders should be drawn. Default  =  TRUE.
  }
  \item{colors}{
     Either "DO", in which case DO colors are supplied or a data.frame with three columns containing the founder letter code, founder strain name and founder color in columns 1:3.
  }
  \item{sex}{
    Character that is either FALSE or M, indicating the sex to use. Only used on X chromosome.
  }
  \item{\dots}{
    Additional arguments to be passed to plot.
  }
}
\value{
  No value is returned.  A plot with the founder coefficients in the top panel and the LOD score in the bottom panel is drawn.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{scanone}}, \code{\link{plot.doqtl}}
}
\examples{
  \dontrun{
    coefplot(qtl, chr = 1)
  }
}
\keyword{ MUGA }
\keyword{ genotyping }
\keyword{ HMM }
