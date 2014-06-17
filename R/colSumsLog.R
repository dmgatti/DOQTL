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
} # colSumsLog
