\name{do.colors}
\alias{do.colors}
\docType{data}
\title{do.colors}
\description{
  The letter codes, strain names and official colors for the DO founders.
}
\usage{do.colors}
\format{
  A data frame with 8 rows and 3 columns.
  \describe{
    \item{\code{CC_Designation}}{a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{FALSE} \code{G} \code{H}}
    \item{\code{Strain}}{a factor with levels \code{129S1/SvImJ} \code{A/J} \code{C57BL/6J} \code{CAST/EiJ} \code{NOD/ShiLtJ} \code{NZO/H1LtJ} \code{PWK/PhJ} \code{WSB/EiJ}}
    \item{\code{R_Color}}{a factor with levels \code{#00A000} \code{#00A0FALSE0} \code{#1010FALSE0} \code{#808080} \code{#9000FALSE4} \code{#FALSE00000} \code{#FALSE08080} \code{#FALSE0FALSE000}}
  }
}
\details{
  This contains the official colors that should be used when plotting data involving the DO founders.
}
\source{
  Copied from \href{http://csbio.unc.edu/CCstatus/index.py?run=AvailableLines.information}{UNC Systems Genetics Website}
}
\examples{
  do.colors
}
\keyword{datasets}
