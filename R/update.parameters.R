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
update.parameters = function(theta, rho, r.t.means, r.t.covars, prsmth,
                    is.founder.F1, founder.means) {

  # Save the dimensions and dimnames of the mean and variance arrays.
  dims = dim(r.t.means)
  namelist = dimnames(r.t.means)

#prsmth = exp(prsmth)
#denom = matrix(0, nrow(r.t.means), dim(r.t.means)[3])
#for(s in 1:nrow(r.t.means)) {
#  r.t.means[,1,] = founder.means$theta
#  denom[s,] = colSums(prsmth[s,,]) + 1
#  numer = colSums(theta * prsmth[s,,])
#  r.t.means[s,1,] = numer / denom[s,]
#  r.t.means[,2,] = founder.means$rho
#  numer = colSums(rho * prsmth[s,,])
#  r.t.means[s,2,] = numer / denom[s,]
#} # for(s)

#for(s in 1:nrow(r.t.means)) {
#  tdiff = theta - matrix(r.t.means[s,1,], nrow(theta), ncol(theta), byrow = T)
#  r.t.covars[s,1,] = colSums(tdiff^2 * prsmth[s,,]) / sqrt(denom[s,])
#  rdiff = rho - matrix(r.t.means[s,2,], nrow(rho), ncol(rho), byrow = T)
#  r.t.covars[s,2,] = colSums(rdiff^2 * prsmth[s,,]) / sqrt(denom[s,])
#} # for(s)


  # Copy the theta & rho variances and covariances into separate matrices.
  tvar = matrix(0, dim(r.t.covars)[1], dim(r.t.covars)[3])
  rvar = matrix(0, dim(r.t.covars)[1], dim(r.t.covars)[3])

  res = .C(update_from_r,
           dims = as.integer(dim(prsmth)),
           t = as.double(theta),
           r = as.double(rho),
           tmeans  = as.double(r.t.means[,1,]),
           rmeans  = as.double(r.t.means[,2,]),
           tvars = as.double(tvar),
           rvars = as.double(rvar),
           prsmth = as.double(prsmth),
           foundertmeans = as.double(founder.means$theta),
           founderrmeans = as.double(founder.means$rho))

  # Reconstruct the data in the format needed by R.
  r.t.means[,1,] = res$tmeans
  r.t.means[,2,] = res$rmeans
  r.t.covars[,1,] = res$tvars
  r.t.covars[,2,] = res$rvars
  
  # Don't let the variances get too small or we will never be able to 
  # place a sample in that state.
  r.t.covars[,1,][r.t.covars[,1,] < 0.001] = 0.001
  r.t.covars[,2,][r.t.covars[,2,] < 0.005] = 0.005

  return(list(r.t.means = r.t.means, r.t.covars = r.t.covars))

} # update.parameters()



