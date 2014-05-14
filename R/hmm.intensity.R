################################################################################
# Generic skeleton for the intensity based HMM.
################################################################################
hmm.intensity = function(data, founders, sex, snps, chr, trans.prob.fxn) {

  # Estimate the (theta, rho) genotype cluster means and variances.
  params = estimate.cluster.params(founders = founders, data = data, chr = chr)
  save(params, file = paste("chr", chr, ".initial.cluster.params.Rata"))

  # Set negative values to NA for this. They will be removed when we take the 
  # means.
  neg = which(founders$theta < 0)
  if(length(neg) > 0) {
    founders$theta[neg] = NA
    founders$rho[neg]   = NA	
  } # if(length(neg) > 0)

  # We need to remove the founders that are not part of the possible 
  # states for the DO/F1 to work.
  founders$theta = founders$theta[founders$code %in% founders$states,]
  founders$rho   = founders$rho[founders$code   %in% founders$states,]
  founders$sex   = founders$sex[founders$code   %in% founders$states]
  founders$code  = founders$code[founders$code  %in% founders$states]

  # Now that we have estimated initial cluster parameters from the founder
  # data, condense the founders down to thier means.
  # For some reason, split and tapply don't work here.
  tmp = by(data = founders$theta, INDICES = founders$code, FUN = colMeans,
           na.rm = TRUE)
  founders$theta = matrix(unlist(tmp), length(founders$states), 
                   ncol(founders$theta), byrow = TRUE, dimnames = 
                   list(founders$states, colnames(founders$theta)))
  tmp = by(data = founders$rho, INDICES = founders$code, FUN = colMeans,
           na.rm = TRUE)
  founders$rho = matrix(unlist(tmp), length(founders$states),  
                 ncol(founders$rho), byrow = TRUE, dimnames = 
                 list(founders$states, colnames(founders$rho)))
  rm(tmp)
  
  # Set NaN values to -1 in the data.
  data$theta[is.na(data$theta) | is.nan(data$theta)] = -1.0
  data$rho[is.na(data$rho) | is.nan(data$rho)] = -1.0
  founders$theta[is.na(founders$theta) | is.nan(founders$theta)] = -1.0
  founders$rho[is.na(founders$rho) | is.nan(founders$rho)] = -1.0

  # Set low theta values to -1 because these tend to represent clusters
  # with low X and Y values that get spread out in the non-Euclidean rho/theta space.
  set.to.neg = which(data$theta < 0.05)
  data$theta[set.to.neg] = -1.0
  data$rho[set.to.neg]   = -1.0
  maxIter = 100
  epsilon = 5 # NOTE: We change this in the first iteration below.
  p = 1
  lastLogLik = -Inf
  logLik = -.Machine$double.xmax
  init.hmm = initialize.hmm(snps = colnames(data$rho), samples = 
             rownames(data$rho), states = founders$states)

  # Estimate the initial emission probabilities.
  b = emission.probs.intensity(data = data, params = params)

  if(attr(data, "sampletype") %in% c("DO", "DOF1")) {
    # For DO samples, we have different transition probabilities for each
    # generation.
    a = trans.prob.fxn(states = founders$states, snps = snps, chr = chr,
        sex = sex, do.gen = data$gen)
  } else {
    # For all other sample types, the transition probabilities are the same
    # for all samples.
    a = list(trans.prob.fxn(states = founders$states, snps = snps, chr = chr,
        sex = sex))
  } # else

  while(p <= maxIter & logLik - lastLogLik > epsilon) {
    print(date())
    print(p)
    # Optional plotting code to watch the HMM progress.
    if(FALSE) {
      layout(matrix(1:2, 1, 2))
      prsmth.plot(1, founders$states, init.hmm$prsmth)
      intensity.mean.covar.plot(s = 3, states = founders$states, theta = data$theta, 
                      rho = data$rho, r.t.means = params$r.t.means, 
                      r.t.covars = params$r.t.covars)
    } # if(plot)
    # Initialize the log-likelihood to a large negative number.
    lastLogLik = logLik
    # Filter and smooth each generation separately because they each have a
    # different transition probability matrix.
    print("Filtering & smoothing...")
    lltmp = matrix(-.Machine$double.xmax, 1, 1)
    for(i in 1:length(a)) {
      # Filter
      gen = 1:nrow(data$theta)
      if(any(names(data) == "gen")) {
        print(paste("gen", names(a)[i]))
        gen = which(data$gen == names(a)[i])
      } # if(any(names(data) == "gen"))
      
      res = .C(C_filter_smooth_intensity, 
               dims = as.integer(dim(init.hmm$prsmth[,gen,,drop = FALSE])),
               a = as.double(a[[i]]), 
               b = as.double(b[,gen,,drop = FALSE]),
               prsmth = as.double(init.hmm$prsmth[,gen,,drop = FALSE]),
               init = as.double(init.hmm$init),
               loglik = as.double(lltmp))
      init.hmm$prsmth[,gen,] = res$prsmth
      lltmp = addLog(lltmp, res$loglik)
    } # for(i)
    logLik = lltmp
    print(paste("LogLik =", logLik))
    # Set the epsilon to be 1/1000 of the starting value.
    # TBD: Is there a better way to set epsilon?
    if(p == 1) {
      epsilon = abs(logLik * 0.001)
    } # if(p == 1)
    # Update the parameters and state means and variances.
    print("update")
    params = parameter.update.intensity(data = data, params = params, 
             prsmth = init.hmm$prsmth, founder.means = founders)
    b = emission.probs.intensity(data = data, params = params)
    p = p + 1
	
    gc()
  } # while(p < maxIter & logLik - lastLogLik > epsilon
  print(date())
  
  # If we reached the end of the iterations.
  if(p >= maxIter) {
    print("Maximum iterations reached")
    r.t.means  = NULL
    r.t.covars = NULL
  } # if(p >= maxIter)
  
  gc()  
  
  return(list(prsmth = init.hmm$prsmth, params = params))
} # hmm.intensity()
