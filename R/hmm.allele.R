################################################################################
# Generic skeleton for the alelle call based HMM.
################################################################################
hmm.allele = function(data, founders, sex, snps, chr, trans.prob.fxn) {

  tmp = em.allele(data = data, founders = founders, sex = sex, snps = snps,
        chr = chr, trans.prob.fxn = trans.prob.fxn)

  return(tmp)

} # hmm.allele

