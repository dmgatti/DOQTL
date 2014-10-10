% -*- coding: utf-8 -*-
\name{NEWS}
\title{News for Package \pkg{DOQTL}}
\encoding{UTF-8}

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
