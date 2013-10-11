################################################################################
# Given the {ACGT} genotypes for a set of samples, create a kinship matrix 
# based on allele sharing.
# Daniel Gatti
# Dan.Gatti@jax.org
# Oct. 27, 2012
################################################################################
# Create a kinship matrix based on allele sharing.
# Arguments: geno: character matrix containing allele calls for each sample and
#                  marker. num.snps x num.samples
kinship.alleles = function(geno) {

  K = matrix(0, ncol(geno), ncol(geno), dimnames = list(colnames(geno),
      colnames(geno)))
  for(i in 1:ncol(geno)) {
    K[,i] = colMeans(geno[,i] == geno)
  } # for(i)

  return(K)

} # kinship.alleles()


# This function takes the 8 state founder probabilities and calculates
# the cosine of the angle between each sample.
# Arguments: probs: 3D array containing founder haplotype contributions.
#                   num.samples x num.founders x num.snps.
#            snps: data.frame containing the locations of the markers.
#                  SNP IDs, chromosomes, Mb & cM locations in columns
#                  1:4. Optional.
#            bychr: boolean that indicates if a separate kinship 
#                   matrix should be made for each chromosome.
kinship.probs = function(probs, snps, bychr = F) {

  K = NULL

  if(bychr) {

    # Make a list of kinship matrices, one for each chromosome, in 
    # which the kinship on one chromosome uses the markers on the 
    # remaining chromosomes.
    if(missing(snps)) {
      stop(paste("snps cannot be missing if bychr is T. Please supply",
           "the marker locations in the snps argument."))
    } # if(missing(snps))

    snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
    probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]
    probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
    snps[,2] = as.character(snps[,2])

    chr = as.character(unique(snps[,2]))
    Kbychr = as.list(chr)
    names(Kbychr) = chr
    keep = as.list(1:length(chr))

    for(i in 1:length(chr)) {
      Kbychr[[i]] = matrix(0, dim(probs)[1], dim(probs)[1], dimnames =
               list(dimnames(probs)[[1]], dimnames(probs)[[1]]))
   
      keep[[i]] = which(snps[,2] == chr[i])
      res = .C(kinship_from_r, 
               dims = as.integer(dim(probs[,,keep[[i]]])),
               probs = as.double(probs[,,keep[[i]]]),
               K = as.double(Kbychr[[i]]))
 
      # Restore the dimensions and dimnames. We have to multiply
      # by the number of SNPs so that the means work out below.
      Kbychr[[i]] = matrix(res$K * length(keep[[i]]), dim(probs)[1],
                    dim(probs)[1], dimnames = 
                    list(dimnames(probs)[[1]], dimnames(probs)[[1]]))
    } # for(i)

    K = as.list(chr)
    names(K) = chr
    num.snps = sapply(keep, length)

    for(i in 1:length(chr)) {
      K[[i]] = matrix(res$K, dim(probs)[1], dim(probs)[1], dimnames = 
               list(dimnames(probs)[[1]], dimnames(probs)[[1]]))
      for(j in chr[chr[i] != chr]) {
        K[[i]] = K[[i]] + Kbychr[[j]]
      } # for(j)
      # Divide by the number of SNPs on the other (!= i) chromosomes.
      K[[i]] = K[[i]] / sum(num.snps[chr != chr[i]])
    } # for(i)

  } else {

    # Use all chromosomes to make a single kinship matrix.
    K = matrix(0, dim(probs)[1], dim(probs)[1], dimnames = list(
        dimnames(probs)[[1]], dimnames(probs)[[1]]))

    res = .C(kinship_from_r, 
             dims = as.integer(dim(probs)),
             probs = as.double(probs),
             K = as.double(K))
 
    # Restore the dimensions and dimnames.
    K = matrix(res$K, dim(probs)[1], dim(probs)[1], dimnames = list(dimnames(probs)[[1]], 
        dimnames(probs)[[1]]))

  } # else

  return(K)

} # kinship.probs()

