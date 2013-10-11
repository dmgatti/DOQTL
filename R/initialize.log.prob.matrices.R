################################################################################
# Initialize the prfilt and prsmthooth matrices.
# Arguments: snps: Vector of SNP IDs.
#            samples: Vector of sample IDs.
#            states: Vector of HMM states.
#            b: the log probability emission matrix.
initialize.log.prob.matrices = function(snps, samples, states, b) {
  
  prpred = array(-.Machine$double.xmax, c(length(states), length(samples),
                 length(snps)), dimnames = list(states, samples, snps))
  prfilt = array(-.Machine$double.xmax, dim(prpred), dimnames =
                 dimnames(prpred))
  prsmth = array(-.Machine$double.xmax, dim(prpred), dimnames =
                 dimnames(prpred))
  
  # Initialize the first SNP of the prpred matrix.
  spl = matrix(unlist(strsplit(states, split = "")), nrow = 2)
  founders = sort(unique(as.vector(spl)))
  homo = which(spl[1,] == spl[2,])
  het  = which(spl[1,] != spl[2,])
  # If there are no hets, then we are probably on the male X Chr.
  # Adjust the probabilities to 1/64 or 1/32.
  if(length(het) == 0) {
    prpred[homo,,1] = log(1.0 / length(founders))
  } else {
    prpred[homo,,1] = log(1.0 / length(founders)^2)
    prpred[het,,1]  = log(2.0 / length(founders)^2)
  } # else

  # Set the row for the founders & F1s that corresponds to the known state
  # in each probability matrix = 1 (log(1) = 0.  Set the rest = 0 
  # (log(0) = -Inf ~= -.Machine$double.xmax)
  # We won't be filtering or smoothing them.
  for(m in 1:dim(prpred)[1]) {
    sam = which(dimnames(prpred)[[1]][m] == dimnames(prpred)[[2]])
    if(length(sam) > 0) {
      # Set the whole sample = 0.
      prpred[,sam,] = -.Machine$double.xmax
      prfilt[,sam,] = -.Machine$double.xmax
      prsmth[,sam,] = -.Machine$double.xmax
      # We know the founder state, so set prob == 1. (log(1) = 0)
      prpred[m,sam,] = 0
      prfilt[m,sam,] = 0
      prsmth[m,sam,] = 0
    } # if(length(sam) > 0)
  } # for(m)
  
  return(list(prpred = prpred, prfilt = prfilt, prsmth = prsmth))

} # initialize.log.prob.matrices()
