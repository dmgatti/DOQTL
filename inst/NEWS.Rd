% -*- coding: utf-8 -*-
\name{NEWS}
\title{News for Package \pkg{DOQTL}}
\encoding{UTF-8}

\section{CHANGES IN VERSION 1.5.3}{
  \subsection{CHANGES}{
    \itemize{
      \item Added get.pgw() to allow users to get genome-wide p-values for LOD or -log10(p-values).
    }
  }
}

\section{CHANGES IN VERSION 1.5.0}{
  \subsection{CHANGES}{
    \itemize{
      \item Major changes to the signficance thresholds. We now produce 
      separate autosomal and X chromosome thresholds. This means that the 
      return value of scanone.perm() is a 2 column matrix (or array).
      The columns contain the autosomal and X Chr maximum statistics, respoectively.
      There is a new function called get.sig.thr() that implements the
      method of Broman et.al., Genetics, 2006 for calculating X chromosome
      thresholds.
    }
  }
}

\section{CHANGES IN VERSION 1.1.8}{
  \subsection{CHANGES}{
    \itemize{
      \item Fixed bug in X chromosome mapping. We now requrire sex as an additive covariate always.
    }
  }
}

\section{CHANGES IN VERSION 1.1.6}{
  \subsection{CHANGES}{
    \itemize{
      \item Moved MUGAExampleData and doMPI to "suggests" field in the DESCRIPTION file.
      \item Changed haplotype reconstruction to work in Euclidean (X,Y) space.
      \item Replaced mclust() with pamk() in initial genotype clustering at each marker. 
            Pamk is faster and the EM algorithm seems to produce equivalent haplotype
            reconstructions.
    }
  }
}



\section{CHANGES IN VERSION 1.1.3}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added support for HS rat haplotype reconstruction.
      \item Added support parallel whole genome association mapping.
    }
  }
}

\section{CHANGES IN VERSION 1.07}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added \code{scanone.assoc} to perform genome wide association mapping.
    }
  }
  \subsection{CHANGES}{
    \itemize{
      \item Fixed bug in \code{assoc.map} where the SDPs were not in the correct order.
    }
  }
}


\section{CHANGES IN VERSION 1.07}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added \code{calc.genoprob2} to run forward/backward algorithm once cluster parameters have been estimated.
    }
  }
  \subsection{CHANGES}{
    \itemize{
      \item Change \code{assoc.map} to use the Sanger VCF file rather than our own custom SNP file.
    }
  }
}



\section{CHANGES IN VERSION 1.06}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added limited support for association mapping in Heterogeneous Stock mice.
    }
  }
  \subsection{CHANGES}{
    \itemize{
      \item \code{calc.genoprob.intensity} Fixed bug in which probs were not written out on the X chromosome when mapping with males only.
    }
  }
}

\section{CHANGES IN VERSION 0.99.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added limited support for Heterogeneous Stock mice.
	  
	  \item \code{assoc.plot} Added strain distribution patterns above the association mapping plot.
    }
  }
  \subsection{CHANGES}{
    \itemize{
      \item \code{rankZ} Fixed bug relating to NA values.
    }
  }
}


\section{CHANGES IN VERSION 0.99.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{read.vcf} reads Sanger SNP VCF files.
	  
	  \item \code{assoc.map} imputes Sanger SNPs onto DO genomes and performs association mapping.
	  
	  \item Fixed bug in \code{kinship.probs} in which kinship per chromosome was not calculated correctly.
	  
	  \item Improved gene layout algorithm in \code{gene.plot}.
    }
  }  
  
  \subsection{CHANGES}{
    \itemize{
      \item \code{scanone} returns p-values and -log10(p-values).
	  
	  \item \code{doqtl.plot} plots either LOD or -log10(p-values).
    }
  }
}
