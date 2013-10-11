estimate.state.locations <-
function(theta, rho, theta.rho.means, theta.rho.covars, nan.snps) {
  warning(paste("Genotype state means not estimated for all states at SNPs", 
          paste(dimnames(theta.rho.means)[[2]][nan.snps], collapse = " "), 
          sep  = " "))

  retval = matrix(0, length(nan.snps), 4, dimnames =
             list(NULL, c("theta.mean", "rho.mean", "theta.var", "rho.var")))
  
  for(s in nan.snps) {
    # Determine which states are missing.
    missing = which(is.nan(theta.rho.means[,s,1]))
    # Get the parental locations.
    spl = strsplit(names(missing), split = "")
    for(i in 1:length(spl)) {
      spl[[i]] = paste(spl[[i]], spl[[i]], sep = "")
      theta.means = theta.rho.means[spl[[i]],s,1]
      rho.means   = theta.rho.means[spl[[i]],s,2]
      retval[i,1] = mean(theta.means)
      retval[i,2] = mean(rho.means)
    } # for(i)
  } # for(s)

  # Simply use the overall mean of the variances.
  retval[,3] = mean(theta.rho.covars[,,1], na.rm = T)
  retval[,4] = mean(theta.rho.covars[,,2], na.rm = T)
  return(retval)
}