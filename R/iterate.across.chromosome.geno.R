################################################################################
# Filter and Smooth back and forth across the current Chr.
# We have to break up the filtering and smoothing by DO generation because
# the transition probabilities are different at each generation.
iterate.across.chromosome.geno = function(geno, a, b, prpred, prfilt,
                                 prsmth, states, plot = F, do.gen) {

  print(date())

  # Filter and smooth each generation separately because they each have a
  # different transition probability matrix.
  for(i in 1:length(a)) {

    # Filter
    print("filter")
    gen = which(do.gen == names(a)[i])
    tmp = filter(a[[i]], b[,gen,], prpred[,gen,], prfilt[,gen,], 
                 is.founder.F1 = rep(F, length(gen)))
    prpred[,gen,] = tmp$prpred
    prfilt[,gen,] = tmp$prfilt
    rm(tmp)

    # Smooth
    print("smooth")
    prsmth[,gen,] = smooth(a[[i]], prpred[,gen,], prfilt[,gen,], prsmth[,gen,],
                             is.founder.F1 = rep(F, length(gen)))

    # Optional plotting code to watch the HMM progress.
    if(plot) {
      plot.prsmth(1, states, prsmth)
    } # if(plot)

  } # for(i)

  print(date())
 
  return(prsmth = prsmth)

} # iterate.across.chromosome.geno()

