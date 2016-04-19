\name{interpolate.markers}
\alias{interpolate.markers}
\title{
  interpolate haplotype or genotype probabilities from one set of markers to another.
}
\description{
  Given a matrix of data with markers in rows and data in columns, interpolate the values in each column from the marker spacing in 'from' to the marker spacing in 'to'.
}
\usage{
  interpolate.markers(data, from, to)
}
\arguments{
  \item{data}{
    Numeric matrix with markers in rows and data in columns. Marker names must be in rownames and must match column 1 of from.
  }
  item{from}{
    Data frame containing at least three columns with marker names, chromosome and position in columns 1, 2 & 3, respectively.
  }
  item{to}{
    Data frame containing at least three columns with marker names, chromosome and position in columns 1, 2 & 3, respectively.
  }
}
\details{
  This function assumes an even and overlapping set of markers in from and to. It may not behave correctly if the marker sets do not overlap.
}
\value{
  Numeric matrix of values interpolated onto the spacing provided in 'to'.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{
    interpolate.markers(data, from, to)
  }
}
