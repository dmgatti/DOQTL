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
emission.probs.intensity2 = function(data, params) {
  # Create a large vector to pass down to C.
  retval = array(0, c(dim(params$x.y.means)[1], nrow(data$x),
                 dim(params$x.y.means)[3]), dimnames = 
                 list(dimnames(params$x.y.means)[[1]], rownames(data$x),
                 dimnames(params$x.y.means)[[3]]))
  
  # Set missing data to have a very low value.
  x = data$x
  x[is.na(data$x)] = -100.0
  y = data$y
  y[is.na(data$y)] = -100.0
  res = .C(C_emission_prob2,
           dims = as.integer(dim(retval)),
           x = as.double(x),
           y = as.double(y),
           xmeans = as.double(params$x.y.means[,1,]),
           ymeans = as.double(params$x.y.means[,2,]),
           xvars  = as.double(params$x.y.covars[,1,]),
           yvars  = as.double(params$x.y.covars[,2,]),
           covars = as.double(params$x.y.covars[,3,]),
           probs  = as.double(retval))
  # Set -Inf values to a the minimum double value.
  res$probs[is.infinite(res$probs)] = -.Machine$double.xmin
  # Create the 3D array of emission probabilities.
  retval = array(res$prob, c(dim(params$x.y.means)[1], nrow(data$x),
           dim(params$x.y.means)[3]), dimnames = 
           list(dimnames(params$x.y.means)[[1]], rownames(data$x),
           dimnames(params$x.y.means)[[3]]))
  return(retval)
} # emission.probs.intensity2()
