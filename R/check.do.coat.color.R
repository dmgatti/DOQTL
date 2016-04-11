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
  mk = which.min(abs(markers[,3] - mean(87.427405, 87.493411)))
  albino = check.genotypes(probs[,,mk[1,1]], coat == "albino", 
           contrib.founders = c("A", "D"))

  # Black: Chr2: 154791402 - 155051012
  mk = markers[markers[,2] == 2,]
  mk = which.min(abs(markers[,3] - mean(154.791402, 155.051012)))
  black = check.genotypes(probs[,,mk[1,1]], coat == "black", 
          contrib.founders = c("A", "B"))

  return(cbind(albino, black))

} # check.do.coat.color()


# Args: probs: Either an 8 or 36 state probs matrix.
#       coat: Booean vector in which TRUE indicates that the mouse has the coat
#             color phenotype. i.e. albino = TRUE
#       contrib.founders: Character vector containing the letter codes for
#             the founders that contribute the coat color allele. i.e. "A", "D".
check.genotype = function(probs, coat, contrib.founders) {

  # Get the two letter genotype calls at the given marker.
  gt = get.max.geno(probs[,,marker,drop = FALSE])
  
  # Compute the genotype combinations that can create the coat color.
  ok.gt = sort(outer(contrib.founders, contrib.founders, paste0))

  results = cbind(sample = name(gt), genotype = gt, coat = coat, 
            gt.ok, coat.ok)

  # Check that the coat color phenotype has the correct genotype.
  results$gt.ok = coat[gt %in% ok.gt]

  # Check that the coat color genotype had the correct coat color.
  results$coat.ok = gt[coat] %in% ok.gt

  return(results)

} # check.genotype()

