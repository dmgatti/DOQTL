################################################################################
# Initialize the prfilt and prsmooth matrices.
# Arguments: snps: Vector of SNP IDs.
#            samples: Vector of sample IDs.
#            states: Vector of HMM states.
#            b: the log probability emission matrix.
# Returns the initial values for the CC or DO.
initialize.hmm = function(snps, samples, states) {
  prsmth = array(-.Machine$double.xmax, c(length(states), length(samples),
                 length(snps)), dimnames = list(states, samples, snps))
  init = rep(0, length(states))
  if(all(nchar(states) == 1)) {
    spl = matrix(unlist(strsplit(states, split = "")), nrow = 1)
    founders = sort(unique(as.vector(spl)))
    init = rep(1.0 / length(founders), length(founders))
  } else {
    # Create the init values.
    spl = matrix(unlist(strsplit(states, split = "")), nrow = 2)
    founders = sort(unique(as.vector(spl)))
    homo = which(spl[1,] == spl[2,])
    het  = which(spl[1,] != spl[2,])
    # If there are no hets, then we are probably on the male X Chr or
    # in the CC.
    # Adjust the probabilities.
    if(length(het) == 0) {
      init[homo] = log(1.0 / length(founders))
    } else {
      init[homo] = log(1.0 / length(founders)^2) # 1/64 for DO
      init[het]  = log(2.0 / length(founders)^2) # 1/32 for DO
    } # else
  } # else
  names(init) = states
  return(list(prsmth = prsmth, init = init))
} # initialize.hmm()
