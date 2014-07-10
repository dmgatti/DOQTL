\name{html.report}
\alias{html.report}
\title{
  Create an HTML report for a set of QTL
}
\description{
  Given a list of QTL objects, create a set of HTML pages summarizing the QTL.
}
\usage{
  html.report(path, qtl, perms, assoc = FALSE)
}
\arguments{
  \item{path}{
    Character string with path to which to write.
  }
  \item{qtl}{
    List containing DOQTL objects as produced by \code{\link{scanone}}.
  }
  \item{perms}{
    Numeric vector containing permutation results for the phenotypes.
  }
  \item{assoc}{
    Boolean, if TRUE, look for corresponding *.Rdata files containing the names of the qtl in the current working directory and plot the association plots. If FALSE (default), do not plot association mapping analysis.
  }
}
\details{
  This function summarizes a set of QTL scans made using \code{\link{scanone}}. FALSEor each QTL in the qtl list, it creates an HTML page with a QTL plot and coefficient plots for the significant QTL. The permutations are used to assess significance. It then summarizes all of these in a table and creates a summary HTML page and a *.csv file with the QTL results.
}
\value{
  No value is returned. HTML report of QTL is written out to the specified directory.
}
\author{
  Daniel Gatti
}
\examples{
   \dontrun{ html.report(path, qtl, perms, merge = FALSE) }
}
\keyword{ HTML }
\keyword{ QTL }
