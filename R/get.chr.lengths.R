#############################################################################
# Get the mouse chromosome lengths.
get.chr.lengths = function(genome = "mm10") {

    genome = get(paste0("BSgenome.Mmusculus.UCSC.", genome))
    
    chrlen = seqlengths(genome)
    chrlen = chrlen[-grep("random|Un", names(chrlen))]
    names(chrlen) = sub("chr", "", names(chrlen))
    chrlen = chrlen[c(1:19, "X", "Y", "M")]
    chrlen = chrlen * 1e-06
    
    return(chrlen)
  
} # get.chr.lengths()
