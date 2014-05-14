################################################################################
# Given the theta and rho intensities and the theta and rho means and 
# covariances, calculate the emission probability from each cluster.
# Arguments: theta: numeric matrix of theta intensities.
#            rho: numeric matrix of rho intensities.
#            r.t.means: 
#            r.t.covars:
################################################################################
emission.probs.intensity = function(data, params) {
  # Create a large vector to pass down to C.
  retval = array(0, c(dim(params$r.t.means)[1], nrow(data$theta),
                 dim(params$r.t.means)[3]), dimnames = 
                 list(dimnames(params$r.t.means)[[1]], rownames(data$theta),
                 dimnames(params$r.t.means)[[3]]))
  
  # Set missing data to have a very low value.
  theta = data$theta
  theta[is.na(data$theta)] = -100.0
  rho = data$rho
  rho[is.na(data$rho)] = -100.0
  
  res = .C(C_emission_prob,
           dims = as.integer(dim(retval)),
           theta = as.double(theta),
           rho   = as.double(rho),
           thetameans = as.double(params$r.t.means[,1,]),
           rhomeans   = as.double(params$r.t.means[,2,]),
           thetavars = as.double(params$r.t.covars[,1,]),
           rhovars   = as.double(params$r.t.covars[,2,]),
           prob = as.double(retval))
  # Create the 3D array of emission probabilities.
  retval = array(res$prob, c(dim(params$r.t.means)[1], nrow(data$theta),
           dim(params$r.t.means)[3]), dimnames = 
           list(dimnames(params$r.t.means)[[1]], rownames(data$theta),
           dimnames(params$r.t.means)[[3]]))
  return(retval)
} # emission.probs.intensity()
