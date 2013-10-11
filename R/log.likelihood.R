################################################################################
# Calculate the minus log-likelihood of the whole model.
# Dependendies: b, prpred, num.states & num.snps.
# Returns: log-likelihood of the entire model.
log.likelihood.intensity = function(b, prpred) {

  retval = 0.0

  for(s in 1:dim(b)[3]) {
    for(i in 1:dim(b)[2]) {
      retval = retval + sum(exp(b[,i,s]) * exp(prpred[,i,s]))
    }
  }
  retval = -log(retval)

#  ll = matrix(0.0, 1, 1)
#  retval = .C(loglik_intensity,
#             dim = as.integer(dim(b)),
#             b = as.double(b),
#             prpred = as.double(prpred),
#             ll = as.double(ll))

#  return(retval$ll)

  return(retval)

} # log.likelihood.intensity()


################################################################################
# Calculate the minus log-likelihood of the whole model.
# Dependendies: b, geno, prpred.
# Returns: log-likelihood of the entire model.
log.likelihood.alleles = function(b, geno, prpred) {

  ll = matrix(0.0, 1, 1)
  res = .C(loglik_alleles, 
           dims = as.integer(c(dim(prpred), dim(b)[1])),
           b = as.double(b),
           geno = as.integer(geno),
           prpred = as.double(prpred),
           ll = as.double(ll))
    
  return(res$ll)

} # log.likelihood.alleles()

