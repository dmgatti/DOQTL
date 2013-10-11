################################################################################
# Generic skeletion for the intensity based HMM.
################################################################################
hmm.intensity = function(data, founders, sex, snps, chr, trans.prob.fxn) {

  # Estimate the (theta, rho) genotype cluster means and variances.
  params = estimate.cluster.params(founders = founders, data = data, chr = chr)

  # Set negative values to NA for this. They will be removed when we take the 
  # means.
  neg = which(founders$theta < 0)
  if(length(neg) > 0) {
    founders$theta[neg] = NA
    founders$rho[neg]   = NA	
  } # if(length(neg) > 0)

  # Now that we have estimated initial cluster parameters from the founder
  # data, condense the founders down to thier means.
  # For some reason, split and tapply don't work here.
  tmp = by(data = founders$theta, INDICES = founders$code, FUN = colMeans,
           na.rm = T)
  founders$theta = matrix(unlist(tmp), length(founders$states), 
                   ncol(founders$theta), byrow = T, dimnames = 
                   list(founders$states, colnames(founders$theta)))
  tmp = by(data = founders$rho, INDICES = founders$code, FUN = colMeans,
           na.rm = T)
  founders$rho = matrix(unlist(tmp), length(founders$states),  
                 ncol(founders$rho), byrow = T, dimnames = 
                 list(founders$states, colnames(founders$rho)))
  rm(tmp)
  tmp = em.intensity(data = data, founders = founders, sex = sex, snps = snps,
        chr = chr, trans.prob.fxn = trans.prob.fxn, params = params)
  gc()
		
  return(list(prsmth = tmp$prsmth, params = tmp$params))

} # hmm.intensity()
