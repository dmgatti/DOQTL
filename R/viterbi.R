################################################################################
# Viterbi algorithm for determining genotype probabilities.
# Daniel Gatti
# Dan.Gatti@jax.org
# Nov. 28, 2013
################################################################################
# Calculate genotype probabilities on one chromosome.
# Arguments: data: list containing either x and y intensities or geno, containing
#                allele calls for each sample. Samples in rows and markers in
#                columns. 
#            params: list containing either rho and theta intensity means and 
#                    covars or emission probabilities for allele calls.
#            snps: data.frame containing the SNP ID, chromosome, Mb and cM
#                  locations in columns 1 through 4, respectively.
viterbi = function(data, founders, params, snps) {
  stop("Viteri algorithm not implemented yet.")
  if("geno" %in% names(data)) {
    viterbi.allele(data, founders, params, snps)
  } else {
    viterbi.intensity(data, founders, params, snps)
  } # else
} # viterbi()

# Helper function for allele call HMM.
viterbi.allele = function(data, founders, params, snps) {
  unique.chr = unique(snps[,2])
  for(chr in unique.chr) {
    ss = which(snps[,2] == chr)
    # Calculate the transition probabilities.
    a = do.trans.probs(states = founders$states, snps = snps, chr = chr, 
        sex = data$sex, do.gen = data$gen)
    # Calculate genotype probabilities.
#    res = .C(C_viterbi,
#             dims = ,
#             a = ,
#             v_path = ,
#             v_prob = ))
  } # for(chr)
} # viterbi.allele()


# Helper function forintensity HMM.
viterbi.intensity = function(data, founders, params, snps) {
  unique.chr = unique(snps[,2])
  for(chr in unique.chr) {
    ss = which(snps[,2] == chr)
    # Calculate the transition probabilities.
    a = do.trans.probs(states = founders$states, snps = snps, chr = chr, 
        sex = data$sex, do.gen = data$gen)
    # Calculate genotype probabilities.
#    res = .C(viterbi_from_r,
#             dims = ,
#             a = ,
#             v_path = ,
#             v_prob = ))
  } # for(chr)
} # viterbi.intensity()
