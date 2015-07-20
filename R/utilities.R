################################################################################
# Utility functions for DO Genotype HMM.
# Daniel Gatti
# Dan.Gatti@jax.org
# Jan. 12, 2012
################################################################################
convert.pos.to.bp = function(pos) {
  if(max(pos) <= 200) {
    pos = pos * 1e6
  } # if(max(pos] <= 200))
  return(pos)
}
################################################################################
# Get the number of autosomes.  We assume that any chromosome with a number is
# an autosome.
# Arguments: snps: matrix of SNPs with at least 4 columns. SNP ID, chr, Mb and
#                  cM in columns 1:4.
# Returns: integer: the number of autosomes.
get.num.auto = function(snps) {
  # Turn off the warnings for a moment while we get the number of autosomes.
  old.warn = options("warn")$warn
  options(warn = -1)
  num.auto = max(as.numeric(snps[,2]), na.rm = TRUE)
  options(warn = old.warn)
  return(num.auto)
  
} # get.num.auto()
################################################################################
# Get the DO genotype states.
# Returns: character vector with the 36 sorted DO states.
get.do.states = function() {
  states = outer(LETTERS[1:8], LETTERS[1:8], paste, sep = "")
  states = sort(states[upper.tri(states, diag = TRUE)])
  states
} # get.do.states()


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
