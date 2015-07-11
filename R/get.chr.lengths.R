#############################################################################
# Get the mouse chromosome lengths.
get.chr.lengths = function(genome = "mm10") {

  bsgenome = get(paste0("BSgenome.Mmusculus.UCSC.", genome))

  chrlen = seqlengths(bsgenome)
  chrlen = chrlen[-grep("random|Un", names(chrlen))]
  names(chrlen) = sub("chr", "", names(chrlen))
  chrlen = chrlen[c(1:19, "X", "Y", "M")]
  chrlen = chrlen * 1e-6
  
  return(chrlen)
  
} # get.chr.lengths()
