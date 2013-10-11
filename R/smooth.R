################################################################################
# Reverse pass through the HMM, using the state at SNP t+1 to update the state
# at SNP t. (Eqn. 13 in Churchill)
# Arguments: a: matrix of transition probabilities (num.states x num.states)
#            prpred: array of predictive probabilities (num.samples x num.snps x
#                    num.states)
#            prfilt: array of filtered probabilities (num.samples x num.snps x
#                    num.states)
#            prsmth: array of smoothed probabilities (num.samples x num.snps x
#                    num.states)
# Dependencies: filter() had been run. prsmth, prfilt, a & num.snps.
# Effects: Populates the prsmth array with the smoothed probabilities.
smooth = function(a, prpred, prfilt, prsmth) {

  # Save the dimensions and dimnames before we make the call down to C.
  dimensions = dim(prpred)
  namelist   = dimnames(prpred)

  res = .C(smooth_from_r,
           dims = as.integer(dim(prpred)),
           a = as.double(a),
           prpred = as.double(prpred),
           prfilt = as.double(prfilt),
           prsmth = as.double(prsmth))
  
  # Restore the dimensions & dimnames.
  prsmth = array(res$prsmth, dimensions, dimnames = namelist)

  return(prsmth)

} # smooth()
