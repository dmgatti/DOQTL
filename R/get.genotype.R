################################################################################
# Given the haplotype probabilities and the location of a single SNP,
# use the haplotype probabilities and the SNP to impute the allele calls in 
# DO mice.
# Daniel Gatti
# Dan.Gatti@jax.org
# Dec. 7, 2015
################################################################################
# Arguments: chr: a single character (or number) that is the chromosome.
#            pos: a positive value that is the SNP location in bp of Mb.
#            snp: The strain pattern of the SNP as 0s and 1s. Must be a named
#                 vector with the strain names.
#            markers: The location of the markers in probs.
#            probs: 3D array containing the haplotype probabilities.
get.genotype = function(chr, pos, snp, markers, probs) {

  # See if the position is in bp or Mb.
  if(pos < 200) {
    pos = pos * 1e6
  } # if(pos < 200)

  # Convert the SNP to numbers.
  snp = unlist(snp)
  names(snp) = make.names(sub("_", ".", names(snp)))
  strains = make.names(do.colors[,2])
  snp = snp[strains]
  snp = as.numeric(factor(snp)) - 1

  # Get the slices from the haplotype probs matrix.
  markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]
  markers = markers[markers[,2] == chr,]
  probs = probs[,,markers[,1]]
  markers = markers[max(which(markers[,3] < pos * 1e-6)):min(which(markers[,3] > pos * 1e-6)),]

  # Get the probs for these markers.
  probs = probs[,,markers[,1], drop = FALSE]
  probs = apply(probs, 1:2, mean)

  # Multiply the two matrices and return the result.
  return(probs %*% snp)

} # get.genotype()
