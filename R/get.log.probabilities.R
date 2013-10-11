get.log.probabilities <-
function(theta, rho, theta.rho.means, theta.rho.covars, states) {

  mach.prec = get.machine.precision()

  # This is 1/(2 * pi * sigma.t * sigma.r)
  prefix = log(0.5 / pi / sqrt(theta.rho.covars[,1,]) / 
               sqrt(theta.rho.covars[,2,]))

  retval = array(0, c(length(states), nrow(rho), dim(theta.rho.means)[3]),
                 dimnames = list(states, rownames(rho),
                 dimnames(theta.rho.means)[[3]]))

  # We assume a diagonal covariance matrix, so the density is just the sum of
  # the marginals.
  for(s in 1:dim(retval)[3]) {
    # (t - mu.t)^2 / sigma.t^2
    z.t.sq = outer(theta.rho.means[,1,s], theta[,s], "-")^2 / 
             theta.rho.covars[,1,s]
    # (r - mu.r)^2 / sigma.r^2
    z.r.sq = outer(theta.rho.means[,2,s], rho[,s],   "-")^2 /
             theta.rho.covars[,2,s]
    retval[,,s] = prefix[,s] - (0.5 * (z.t.sq + z.r.sq))

    # Make the sum of each row = 1 by dividing by the sum.
    retval[,,s] = retval[,,s] - rowSumsLog(retval[,,s], mach.prec)
  } # for(s)

  return(retval)

}