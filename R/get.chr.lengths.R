#############################################################################
# Get the mouse chromosome lengths.
get.chr.lengths = function(genome) {

  genome = paste0("BSgenome.Mmusculus.UCSC.", genome)
  obj = get(genome)

  # Get the chromosome lengths.
  chrlen = seqlengths(obj)
  chrlen = chrlen[-grep("Un|random", names(chrlen))]
  chrlen = chrlen * 1e-6

  return(chrlen)

} # get.chr.lengths()
