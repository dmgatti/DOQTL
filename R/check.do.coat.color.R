################################################################################
# Check the DO coat color at the albino and agouti locus.
# Daniel Gatti
# Dan.Gatti@jax.org
# Feb. 2, 2015
# Arguments: probs: Numeric 3D array containing either the diplotype probs
#                   or the haplotype probs. Samples x states x markers.
#            markers: Data.frame containing the marker locations. Marker ID,
#                     chr, Mb and cM postions in columns 1 through 4.
#            coat: Character vector containing one of "agouti", "albino", 
#                  "black" or NA. Sample IDs in names.
# The number and names of samples must match between probs and coat.
################################################################################
check.do.coat.color = function(probs, markers, coat) {

  # Albino: Chr 7: 87427405 - 87493411
  mk = markers[markers[,2] == 7,]
  mk = mk[which.min(abs(mk[,3] - mean(87.427405, 87.493411))),1]
  albino = check.genotype(probs[,,mk], coat == "albino", 
           contrib.founders = c("A", "D"))

  # Black: Chr2: 154791402 - 155051012
  mk = markers[markers[,2] == 2,]
  mk = mk[which.min(abs(mk[,3] - mean(154.791402, 155.051012))),1]
  black = check.genotype(probs[,,mk], coat == "black", 
          contrib.founders = c("A", "B"))

  return(cbind(albino, black))

} # check.do.coat.color()


# Args: pr: Either an 8 or 36 state probs matrix.
#       coat: Booean vector in which TRUE indicates that the mouse has the coat
#             color phenotype. i.e. albino = TRUE
#       contrib.founders: Character vector containing the letter codes for
#             the founders that contribute the coat color allele. i.e. "A", "D".
check.genotype = function(pr, coat, contrib.founders) {

  # Get the two letter genotype calls at the given marker.
  rnk = apply(pr, 1, order)
  gt = cbind(LETTERS[1:8][rnk[7,]], LETTERS[1:8][rnk[8,]])
  gt = t(apply(gt, 1, sort))
  homo = apply(pr, 1, max)
  homo = which(homo > 0.75)
  gt[homo,] = cbind(LETTERS[1:8][rnk[8,homo]], LETTERS[1:8][rnk[8,homo]])
  gt = apply(gt, 1, paste0, collapse = "")
  names(gt) = rownames(pr)
  
  # Compute the genotype combinations that can create the coat color.
  possible.gt = sort(outer(contrib.founders, contrib.founders, paste0))

  results = data.frame(sample = names(gt), genotype = gt, coat = coat, 
            gt.ok = rep(NA, length(gt)), coat.ok = rep(NA, length(gt)))

  # Check that the coat color phenotype has the correct genotype.
  samples.with.color.geno = which(results$genotype %in% possible.gt)
  results$gt.ok[samples.with.color.geno] = results$coat[samples.with.color.geno] == TRUE

  # Check that the coat color genotype had the correct coat color.
  samples.with.coat.color = which(coat)
  results$coat.ok[samples.with.coat.color] = results$genotype[samples.with.coat.color] %in% possible.gt

  return(results)

} # check.genotype()

