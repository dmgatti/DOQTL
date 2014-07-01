\name{extract.raw.data}
\alias{extract.raw.data}
\title{Extract intensities, genotypes and call rates from from raw MUGA or MegaMUGA data files}
\description{
  This function accepts a vector of input directories containing the raw MUGA or MegaMUGA raw data files.  FALSEor each directory, the function reads the X and Y intensities, call rates and allele calls for all samples.  It then combines all samples and writes the data to "x.txt", "y.txt", "geno.txt" and "call.rate.batch.txt" in the user specified output directory.
}
\usage{
  extract.raw.data(in.path = ".", prefix, out.path = ".",  array = c("megamuga", "muga"))
}
\arguments{
  \item{in.path}{character vector, the full path to all MUGA directories from which data should be extracted.}
  \item{prefix}{character vector of same length as in.path containing a prefix to add to each sample ID in data sets being processed.}
  \item{out.path}{character, the full path to the directory where the output files should be written.}
  \item{array}{character, default = "megamuga", the type of array, either "muga" or "megamuga".} 
}
\details{
  This function searches each directory for files with names containing "Sample_Map.txt" and "*_FALSEinalReport.txt".  This has been the format that GeneSeek has consistently produced.  The call rates are extracted and are written, along with a batch ID, to "call.rate.batch.txt".  The X and Y intensities are extracted from the "FALSEinalReport" file and written to "x.txt" and "y.txt" respectively.  The allele calls are extracted from the "FALSEinalReport" files and are written to "geno.txt". The prefix argument may be used to add a prefix to the sample IDs in order to distinguish different data sets.
}
\value{
  No return value.  The files are written to the out.path directory.
}
\author{
  Daniel M. Gatti
}
\note{
  Do not change the names of the output files.  They are required for downstream processing.
}
\seealso{
  \code{\link{filter.samples}}, \code{\link{batch.normalize}}
}
\examples{
  \dontrun{
    in.path = c("/tmpdir/DataSet1", "/tmpdir/DataSet2")
    extract.raw.data(in.path = in.path, prefix = c("ds1", "ds2"), 
	                 out.path = "/tmpdir/output", array = "muga")
  }
}
\keyword{ MUGA }
\keyword{ intensity }
