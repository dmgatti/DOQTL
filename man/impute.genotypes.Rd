\name{impute.genotypes}
\alias{impute.genotypes}
\title{
  Impute Sanger SNPs onto mouse genomes.
}
\description{
  Given a GRanges object containing a single range on one chromosome and haplotype probabilities, impute the Sanger SNPs onto the DO genomes.
}
\usage{
  impute.genotypes(gr, probs, markers, vcf.file, hq = TRUE, 
  cross = c("DO", "CC", "DOF1", "HS", "HSrat", "other"))
}
\arguments{
  \item{gr}{
     object that contains the genomic range in which to impute SNPs. For now, this must be a single, continuous range.
  }
  \item{probs}{
    3D numeric array containing the haplotype probabilities. samples x 8 founders x markers. All dimensions must be contain dimnames.
  }
  \item{markers}{
    data.frame containing at least 3 columns that include the marker ID, chr and postion of each marker. nrow must be the same as dim(probs)[3] and markers[,1] must equal dimnames(probs)[[3]].
  }
  \item{vcf.file}{
    String containing the full path to the Sanger SNP VCF file.
  }
  \item{hq}{
    Boolean indicating whether to use only high quality SNPs. Default = TRUE.
  }
  \item{cross}{
    Character string that is the cross type. One of "DO", "CC", "DOF1", "HS", "HSrat", "other")
  }
}
\details{
  This function takes the mean of the haplotype probabilities between two markers and does not linearly interpolate between markers.
}
\value{
  A file written out to the outfile name.
}
\author{
  Daniel Gatti
}
\examples{
  \dontrun{
    impute.genotypes(gr, probs, markers, vcf.file, hq = TRUE, 
    cross = "DO")
  }
}
