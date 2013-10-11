################################################################################
# Using the smoothed probabilities, update the state means and covariances.
# Arguments: geno: matrix of alelels calls coded as 0 = AA, 1 = het,
#                  2 = BB, 3 = no call. num_samples x num_snps.
#            b: 3D array of emission probabilities. num_symbols x num_states
#               x num_snps.
#            pseudocounts: 3D array of pseudocounts from the founders.
#                          num_symbols x num_states x num_snps.
#            prsmth: 3D array with smoothed state probabilities for each sample.
# Dependencies: filter() and smooth() have been run.  prsmth, b, geno.
# (Eqn. 12 in Churchill)
update.parameters.alleles = function(geno, b, pseudocounts, prsmth) {

  # Save the dimensions and dimnames of the b array.
  dims = dim(b)
  namelist = dimnames(b)

  res = .C(update_alleles_from_r,
           dims = as.integer(c(dim(prsmth), dim(b)[1])),
           geno = as.integer(geno),
           b = as.double(b),
           pseudo = as.double(pseudocounts),
           prsmth = as.double(prsmth))

  # Reconstruct the data in the format needed by R.
  b = array(data = res$b, dim = dims, dimnames = namelist)
  
  return(b)

} # update.parameters.alleles()



