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
parameter.update.alleles = function(geno, b, pseudocounts, prsmth) {
  # Save the dimensions and dimnames of the b array.
  dims = dim(b)
  namelist = dimnames(b)
  res = .C(C_update_alleles,
           dims = as.integer(c(dim(prsmth), dim(b)[1])),
           geno = as.integer(geno),
           b = as.double(b),
           pseudo = as.double(pseudocounts),
           prsmth = as.double(prsmth))
  # Reconstruct the data in the format needed by R.
  b = array(data = res$b, dim = dims, dimnames = namelist)
  
  return(b)
} # parameter.update.alleles()

################################################################################
# Using the smoothed probabilities, update the state means and covariances.
# Arguments: x: matrix of X intensities. num_samples x num_snps.
#            y: matrix of Y intensities. num_samples x num_snps.
#            r.t.means: 3D array of X & Y state intensity means.
#            r.t.covars: 3D array of X & Y state intensity variances.
#            prsmth: 3D array with smoothed state probabilities for each sample.
#            is.founder.F1: boolean vector that is T at indices where the
#                           sample is a founder or F1.
#            founder.means: list, with two matrices named theta and rho.
#                           Each matrix is the mean founder value for each SNP.
# Dependencies: filter() and smooth() have been run.  prsmth, prfilt, b.
# (Eqn. 12 in Churchill)
parameter.update.intensity = function(data, params, prsmth, founder.means) {
  # Save the dimensions and dimnames of the mean and variance arrays.
  dims     = dim(params$r.t.means)
  namelist = dimnames(params$r.t.means)
  res = .C(C_update_intensity,
           dims = as.integer(dim(prsmth)),
           t = as.double(data$theta),
           r = as.double(data$rho),
           tmeans = as.double(params$r.t.means[,1,]),
           rmeans = as.double(params$r.t.means[,2,]),
           tvars  = as.double(params$r.t.covars[,1,]),
           rvars  = as.double(params$r.t.covars[,2,]),
           prsmth = as.double(prsmth),
           foundertmeans = as.double(founder.means$theta),
           founderrmeans = as.double(founder.means$rho))
  # Reconstruct the data in the format needed by R.
  params$r.t.means[,1,]  = res$tmeans
  params$r.t.means[,2,]  = res$rmeans
  params$r.t.covars[,1,] = res$tvars
  params$r.t.covars[,2,] = res$rvars
  
  # Don't let the variances get too small or we will never be able to 
  # place a sample in that state
  params$r.t.covars[,1,][params$r.t.covars[,1,] < 0.001] = 0.001
  params$r.t.covars[,2,][params$r.t.covars[,2,] < 0.005] = 0.005
  return(params)
} # parameter.update.intensity()
