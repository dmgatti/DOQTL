################################################################################
# Functions to add logs on the untransformed scale.  This is not just adding
# logs (as in multiplying the untransformed numbers).  These are used when we
# must operate on the log scale to avoid underflows, but need to sum the values.
# Daniel Gatti
# Dan.Gatti@jax.org
# April 9, 2013
################################################################################
addLog = function(x, y) {
  retval = 0.0
  if(is.infinite(x) && x < 0) {
    retval = y
  } else if(is.infinite(y) && y < 0) { 
    retval =  x
  } else if (x >= y) {
    retval = x + log1p(exp(y - x))
  } else {
    retval = y + log1p(exp(x - y))
  } # else
  return(retval)
} # addLog()
# Add the values in a vector of logs.
addLogVector = function(x) {
  # Get the largest value and save it.
  max.index = which.max(x)
  retval = x[max.index]
  # Divide all values by the largest value.
  x = x - retval
  # Set those values less than the machine precision that 
  # are too small to influence the largest value to NA.
  x[x < get.machine.precision()] = NA
  # Set the largest value to NA because we already have it in retval.
  x[max.index] = NA
  # Return the largest value * the sum of the remaining values.
  return(retval + log1p(sum(exp(x), na.rm = TRUE)))
} # addLogVector()
################################################################################
# Sum the rows of a matrix that is on a log(e) scale, implementing
# log(sum(exp(x.i))) = x.0 + log(sum(exp(x.i - x.0)))
# where the x's are sorted and x.0 is the largest value.
# We remove values that are more than 35 logs less than the maximum value.
# Note that exp(-35) ~ 10^-16.
colSumsLog = function(logmat) {
 
  max.index = apply(logmat, 2, which.max) + 0:(ncol(logmat) - 1) * nrow(logmat)
  retval = logmat[max.index]
  logmat[max.index] = NA
  logmat = logmat - matrix(retval, nrow(logmat), ncol(logmat), byrow = TRUE)
  logmat[logmat < get.machine.precision()] = NA
  retval = retval + log1p(colSums(exp(logmat), na.rm = TRUE))
  return(retval)
} # colSumsLog()
################################################################################
# Sum the rows of a matrix that is on a log(e) scale, implementing
# log(sum(exp(x.i))) = x.0 + log(sum(exp(x.i - x.0)))
# where the x's are sorted and x.0 is the largest value.
# We remove values that are more than 35 logs less than the maximum value.
# Note that exp(-35) ~ 10^-16.
rowSumsLog = function(logmat) {
  max.index = apply(logmat, 1, which.max) + 0:(nrow(logmat) - 1) * ncol(logmat)
  retval = logmat[max.index]
  logmat[max.index] = NA
  logmat = logmat - retval
  logmat[logmat < get.machine.precision()] = NA
  retval = retval + log1p(rowSums(exp(logmat), na.rm = TRUE))
  return(retval)
} # rowSumsLog()
