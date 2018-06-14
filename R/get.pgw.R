################################################################################
# Given a set of permutations and a LOD or -log10(p-value), caluclate the
# genome wide p-value.
# Daniel Gatti
# dan.gatti@jax.org
# Nov. 5, 2015
################################################################################
get.pgw = function(stat, chr, perms) {
  if(length(stat) != length(chr)) {
    stop("The length of the statistics must be the same as the length of the",
         "chromosomes.")
  } # if(length(stat) != length(chr))
  if(min(stat) >= 0 & max(stat) <= 1) {
    warning("The values in stat are all between 0 and 1. If they are p-values,
            please convert them to -log10(p-values).")
  } # if(min(stat) >= 0 & max(stat) <= 1)
  if(min(perms) >= 0 & max(perms) <= 1) {
    warning("The values in perms are all between 0 and 1. If they are p-values,
            please convert them to -log10(p-values).")
  } # if(min(perms) >= 0 & max(perms) <= 1)
  # Get the chromosome lengths.
  chrlen = get.chr.lengths()
  len.auto = sum(chrlen[1:19])
  len.X = chrlen["X"]
  len.all = len.auto + len.X
  # Set up the return value.
  pgw = rep(1, length(stat))
  names(pgw) = names(stat)
  isA = which(chr != "X")
  isX = which(chr == "X")
  # Autosomes
  if(length(isA) > 0) {
    q = sapply(lapply(stat[isA], "<", perms[,"A"]), mean, na.rm = TRUE)
    pgw[isA] = 1.0 - (1.0 - q)^(len.all / len.auto)  
  } # if(length(!isX) > 0)
  # X chromosome.
  if(length(isX) > 0) {
    q = sapply(lapply(stat[isX], "<", perms[,"X"]), mean, na.rm = TRUE)
    pgw[isX] = 1.0 - (1.0 - q)^(len.all / len.X)  
  } # if(length(isX) > 0)
  return(pgw)
} # get.pgw()
