################################################################################
# Filter and Smooth back and forth across the current Chr.
# We have to break up the filtering and smoothing by DO generation because
# the transition probabilities are different at each generation.
em.intensity = function(data, founders, sex, snps, chr, trans.prob.fxn, 
               params, plot = F) {

  # Set NaN values to -1 in the data.
  data$theta[is.nan(data$theta)] = -1.0
  data$rho[is.nan(data$rho)] = -1.0
  maxIter = 100

  epsilon = 5 # NOTE: We change this in the first iteration below.
  p = 1
  lastLogLik = -Inf
  logLik = -.Machine$double.xmax

  init.hmm = initialize.hmm(snps = colnames(data$rho), samples = 
             rownames(data$rho), states = founders$states)

  # Estimate the initial emission probabilities.
  b = emission.probs.intensity(data = data, params = params)

  if(attr(data, "sampletype") == "DO") {
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
    if(F) {
      layout(matrix(1:2, 1, 2))
      plot.prsmth(1, founders$states, init.hmm$prsmth)
      plot.mean.covar(s = 3, states = founders$states, theta = data$theta, 
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
      
      res = .C(filter_smooth_intensity, 
               dims = as.integer(dim(init.hmm$prsmth[,gen,,drop = F])),
               a = as.double(a[[i]]), 
               b = as.double(b[,gen,,drop = F]),
               prsmth = as.double(init.hmm$prsmth[,gen,,drop = F]),
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
    params = update.parameters.intensity(data = data, params = params, 
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
    
  return(list(params = params, prsmth = init.hmm$prsmth))

} # em.intensity()

