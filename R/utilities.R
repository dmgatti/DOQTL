################################################################################
# Utility functions for DO Genotype HMM.
# Daniel Gatti
# Dan.Gatti@jax.org
# Jan. 12, 2012
################################################################################
convert.pos.to.bp = function(pos) {
  if(max(pos) <= 300) {
    pos = pos * 1e6
  } # if(max(pos] <= 300))
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
  genome.prefix = substring(genome, 1, 2)
  if(genome.prefix != "mm" & genome.prefix != "rn") {
    stop(paste("get.chr.lengths: Unrecognized genome:", genome, ".\nDOQTL",
         "currently support the mouse and rat genomes (i.e. mm10 or rn6)"))
  } # if(genome.prefix != "mm" & genome.prefix != "rn")
  species = "Mmusculus"
  if(genome.prefix == "rn") {
    species = "Rnorvegicus"
  } # if(substring(genome, 1, 2) == "rn")
  genome = get(paste0("BSgenome.", species, ".UCSC.", genome))
  chrlen = seqlengths(genome)
  chrlen = chrlen[-grep("random|Un", names(chrlen))]
  names(chrlen) = sub("chr", "", names(chrlen))
  
  old.warn = options("warn")$warn
  options(warn = -1)  
  autosomes = sort(which(!is.na(as.numeric(names(chrlen)))))
  options(warn = old.warn)
  chrlen = chrlen[c(autosomes, "X", "Y", "M")]
  chrlen = chrlen * 1e-6
  return(chrlen)
  
} # get.chr.lengths()
#############################################################################
# Get the precision on this machine.
get.machine.precision = function() {
  return(log10(.Machine$double.eps) / log10(exp(1)))
}
