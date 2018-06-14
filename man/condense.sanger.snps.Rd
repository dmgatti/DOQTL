\name{condense.sanger.snps}
\alias{condense.sanger.snps}
\title{
  Create an HDF5 file with the unique SNP patterns between each pair of markers.
}
\description{
  Given a set of markers, a SNP VCF and a set of founder strains that are in the Sanger Mouse Genomes, create an HDF5 file that contains one object for each pair of markers.
}
\usage{
  condense.sanger.snps = function(markers, snp.file, strains, hdf.file, ncl = 1)
}
\arguments{
  \item{markers}{
    Data.frame containing the markers locations. There must be at least three columns containing the marker name, the chromosome and the Mb position. Other columns are ignored.
  }
  \item{snp.file}{
    Character string containing the full path to the Sanger SNP VCF file. Get the file from ftp://ftp-mouse.sanger.ac.uk/.
  }
  \item{strains}{
    Character vector containing the names of strains to extract from the Sanger VCF file. Obtain the names using \code{\link{get.vcf.strains}}. 
  }
  \item{hdf.file}{
    Character string containing the name of the HDF5 file. NOTE: this function will not overwrite an existing file. 
  }
  \item{ncl}{
    Integer that is the number of nodes to use for parallel execution. Default = 1. 
  }
}
\details{
  The markers are usually MUGA series markers, but they can be any set of genomic positions. The markers are broken up by chromosome and written out to and HDF5 file in which each data object is named using the proximal and distal marker named separated by and "_". At the beginning of the chromosome, the object is named using "start" and the first marker separated by an "_". At the end of the chromosome, the object is named using the last marker and "end" separated by an "_".
  Each object is a list containing three elements:
  sdps: Numeric matrix containing the unique strain distribution patterns in the founder strains.
  map:  Numeric vector indicating the index of the SNP at which each SDP occurs with the current pair of markers.
  snps: Data.frame containing SNP names, positions, reference and alternate alleles.
  NOTE: We assume (incorrectly) that all SNPs are bimorphic. SNPs are coded as 0 for "equal to reference" and 1 for "not equal to reference".
}
\value{
  Writes out and HDF5 file.
}
\author{
  Daniel Gatti
}
\seealso{
  \code{\link{read.vcf}}, \code{\link{read.vcf}}
}
\examples{
  \dontrun{ 
    load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
	vcf.file = "mgp.v4.snps.dbSNP.vcf.gz"
    condense.sanger.snps(markers = muga_snps, snp.file = vcf.file, strains = get.vcf.strains(vcf.file), hdf.file = "output.h5") 
  }
}
\keyword{ Sanger }
