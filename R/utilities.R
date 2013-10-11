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
# Convert the text based smoothed probability files to Rdata objects.
# Reads in the smoothed probability text files and writes them out as *.Rdata
# objects.
# Arguments: path: character, full path to the directory with the
#                  *.genotype.probs.txt files.
# Returns: Nothing.
convert.geno.prob.to.Rdata = function(path) {

  files = dir(path, pattern = "genotype.probs.txt", full.names = T)
  for(f in files) {
    print(f)
    prsmth = as.matrix(exp(read.delim(f)))
    save(prsmth, file = sub("txt", "Rdata", f))
  } # for(f)
  
} # convert.geno.prob.to.Rdata()


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
  num.auto = max(as.numeric(snps[,2]), na.rm = T)
  options(warn = old.warn)

  return(num.auto)
  
} # get.num.auto()


################################################################################
# Get the DO genotype states.
# Returns: character vector with the 36 sorted DO states.
get.do.states = function() {
  states = outer(LETTERS[1:8], LETTERS[1:8], paste, sep = "")
  states = sort(states[upper.tri(states, diag = T)])
  states
} # get.do.states()
