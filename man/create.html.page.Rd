\name{create.html.page}
\alias{create.html.page}
\title{
  Create an HTML QTL report
}
\description{
  Given a DOQTL object, create an HTML page that reports the QTL and creates QTL plots. Permutations for assessing significance thresholds can be supplied.
}
\usage{
  create.html.page(path, qtl, pheno.name, perms)
}
\arguments{
  \item{path}{
    Character string containing the path to which to write the report.
  }
  \item{qtl}{
    DOQTL object containing a list with two elemends: lod and coef.
  }
  \item{pheno.name}{
    Character string containing the phenotype name.
  }
  \item{perms}{
    Numeric vector containing the permutation LOD scores for this phenotype.
  }
}
\details{
  The function creates an HTML page with a QTL plot, a table of significant QTL and coefficient plots for the significant loci.
}
\value{
  Data.frame with the significant QTL from this phenotype.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{html.report}}
}
\examples{
  \dontrun{ create.html.page(path, qtl, pheno.name, perms) }
}
\keyword{ HTNL }
