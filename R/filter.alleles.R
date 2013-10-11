################################################################################
# Forward pass through the HMM, using the state at SNP t-1 to update the state
# at SNP t. (Eqns. 10 & 11 in Churchill, 1989)
# Arguments: geno: numeric matrix of genotype allele calls. 0 = AA, 1 = het,
#                  2 = BB, 3 = no call.
#            a: matrix of transition probabilities (num.states x num.states)
#            b: array of emission probabilities (num.alleles x num.states x
#               num.samples)
#            prpred: array of predictive probabilities (num.snps x num.states x
#                    num.samples)
#            prfilt: array of filtered probabilities (num.snps x num.states x
#                    num.samples)
# Dependencies: a, b, prpred, prfilt.
# Effects: Populates the prpred and prfilt arrays with the predictive and
# filtered probabilities.
filter.alleles = function(geno, a, b, prpred, prfilt) {

  # Save the dimensions and dimnames before we make the call down to C.
  dimensions = dim(prpred)
  namelist   = dimnames(prpred)

  res = .C(filter_alleles_from_r, 
            dims = as.integer(c(dim(prpred), dim(b)[1])),
            geno = as.integer(geno),
            a = as.double(a), 
            b = as.double(b),
            prpred = as.double(prpred),
            prfilt = as.double(prfilt))

  # Restore the dimensions and dimnames.
  prpred = array(res$prpred, dimensions, dimnames = namelist)
  prfilt = array(res$prfilt, dimensions, dimnames = namelist)

  return(list(prpred = prpred, prfilt = prfilt))

} # filter.alleles()
