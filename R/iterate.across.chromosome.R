################################################################################
# Filter and Smooth back and forth across the current Chr.
# We have to break up the filtering and smoothing by DO generation because
# the transition probabilities are different at each generation.
iterate.across.chromosome = function(theta, rho, a, b, prpred, prfilt, prsmth,
                                     r.t.means, r.t.covars, is.founder.F1,
                                     states, plot = F, do.gen, founder.means) {

  mach.prec = get.machine.precision()
  maxIter = 100
  epsilon = 5 # NOTE: We change this in the first iteration below.
  p = 1
  lastLogLik = -Inf
  logLik = -.Machine$double.xmax

  while(p < maxIter & logLik - lastLogLik > epsilon) {
    print(date())
    print(p)

    # Optional plotting code to watch the HMM progress.
    if(plot) {
      layout(matrix(1:2, 1, 2))
      plot.prsmth(sum(is.founder.F1) + 11, states, prsmth)
      plot.mean.covar(s = 11, states = states, theta = theta, rho = rho,
                      t.r.means = r.t.means, t.r.covars = r.t.covars)
    } # if(plot)

    # Filter and smooth each generation separately because they each have a
    # different transition probability matrix.
    for(i in 1:length(a)) {

      # Filter
      print("filter")
      gen = which(do.gen == names(a)[i])
      tmp = filter.intensity(a[[i]], b[,gen,], prpred[,gen,], prfilt[,gen,], 
            is.founder.F1 = rep(F, length(gen)))
      prpred[,gen,] = tmp$prpred
      prfilt[,gen,] = tmp$prfilt
      rm(tmp)

      # Smooth
      print("smooth")
      prsmth[,gen,] = smooth(a[[i]], prpred[,gen,], prfilt[,gen,],
                      prsmth[,gen,], is.founder.F1 = rep(F, length(gen)))
    } # for(i)

    # Calculate log-likelihood of the model
    lastLogLik = logLik
    logLik = log.likelihood.intensity(b = b, prpred = prpred, 
             is.founder.F1 = is.founder.F1)
    print(logLik)
    
    if(p == 1) {
      epsilon = abs(logLik * 0.001)
    } # if(p == 1)

    # Update the parameters and state means and variances.
    print("update")
    tmp = update.parameters(theta, rho, r.t.means, r.t.covars, prsmth,
          is.founder.F1, founder.means)
    r.t.means  = tmp$r.t.means
    r.t.covars = tmp$r.t.covars
    rm(tmp)

    b = get.emission.probs(theta, rho, r.t.means, r.t.covars)
    p = p + 1
    print(date())

  } # while(p < maxIter & logLik - lastLogLik > epsilon

  # If we reached the end of the interations.
  if(p >= maxIter) {
    print("Maximum iterations reached")
    r.t.means  = NULL
    r.t.covars = NULL
  } # if(p >= maxIter)
    
  return(list(r.t.means = r.t.means, r.t.covars = r.t.covars, prsmth = prsmth))

} # iterate.across.chromosome()

