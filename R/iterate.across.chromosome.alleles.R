################################################################################
# Filter and Smooth back and forth across the current Chr.
# We have to break up the filtering and smoothing by DO generation because
# the transition probabilities are different at each generation.
iterate.across.chromosome.alleles = function(geno, a, b, prpred, prfilt,
                                 prsmth, states, plot = F, do.gen) {

  mach.prec = get.machine.precision()
  maxIter = 100
  epsilon = 5 # NOTE: We change this in the first iteration below.
  p = 1
  lastLogLik = -Inf
  logLik = -.Machine$double.xmax

  # Save the initial emission probabilities as pseudocounts.
  pseudocounts = exp(b)

  while(p < maxIter & logLik - lastLogLik > epsilon) {
    print(date())
    print(p)

    # Optional plotting code to watch the HMM progress.
    if(plot) {
      plot.prsmth(2, states, prsmth)
    } # if(plot)

    # Filter and smooth each generation separately because they each have a
    # different transition probability matrix.
    print("Running HMM ...")
    for(i in 1:length(a)) {

      # Filter
      gen = which(do.gen == names(a)[i])
      tmp = filter.alleles(geno = geno[gen,], a = a[[i]], b = b,
            prpred = prpred[,gen,], prfilt = prfilt[,gen,])
      prpred[,gen,] = tmp$prpred
      prfilt[,gen,] = tmp$prfilt
      rm(tmp)

      # Smooth
      prsmth[,gen,] = smooth(a[[i]], prpred[,gen,], prfilt[,gen,],
                      prsmth[,gen,], is.founder.F1 = rep(F, length(gen)))
    } # for(i)

    # Calculate log-likelihood of the model
    lastLogLik = logLik
    logLik = log.likelihood.alleles(b = b, geno = geno, prpred = prpred, 
             is.founder.F1 = is.founder.F1)
    print(logLik)
    
    if(p == 1) {
      epsilon = abs(logLik * 0.001)
    } # if(p == 1)

    # Update the parameters and state means and variances.
    print("Updating Parameters ...")
    tmp = update.parameters.alleles(geno = geno, b = b, 
          pseudocounts = pseudocounts, prsmth = prsmth)
    b = tmp
    rm(tmp)

    p = p + 1

  } # while(p < maxIter & logLik - lastLogLik > epsilon

  print(date())

  # If we reached the end of the interations.
  if(p >= maxIter) {
    print("Maximum iterations reached")
  } # if(p >= maxIter)
    
  return(prsmth)


} # iterate.across.chromosome.alleles()

