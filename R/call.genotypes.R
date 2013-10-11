################################################################################
# Given the HMM parameters, calculate the Viterbi path through the genotype
# states for each sample.
# Daniel Gatti
# Dan.Gatti@jax.org
# May, 2, 2012
################################################################################
# Arguments: cross: list with 
call.genotypes = function(cross, theta.means, rho.means, theta.rho.covars) {

  data(states)

  # Get the chromosomes.
  data(snps)
  snps = snps[snps[,1] %in% colnames(cross$x),]
  unique.chr = unique(snps[,2])

  # Create rho and theta values.
  rho = cross$x + cross$y
  theta = (2.0/pi) * atan2(cross$y, cross$x)

  # Extract sex and do generation number.
  sex = cross$sex.gen[,1]
  gen = cross$sex.gen[,2]
  rm(cross)

  # Go through each autosome.
  autosomes = unique.chr[unique.chr %in% 1:19]
  retval = list()
#  for(chr in autosomes) {
chr = 1
print(chr)
    cur.snps = which(snps[,2] == chr)

    theta.rho.means = array(0, c(nrow(rho.means), 2, length(cur.snps)),
                      dimnames = list(rownames(rho.means), c("theta", "rho"),
                      colnames(rho.means)[cur.snps]))
    theta.rho.means[,1,] = theta.means[,cur.snps]
    theta.rho.means[,2,] = rho.means[,cur.snps]

    # Calculate the transition probabilities.
    a = create.log.transition.matrices(states, snps[cur.snps,],
        gen, chr = as.character(chr), sex = "F")

    # Calculate the emission probabilities.
    b = get.log.probabilities(theta[,cur.snps], rho[,cur.snps],
        theta.rho.means, theta.rho.covars, states)

    # Calculate each DO generation separately because each one has different
    # transition matrices.
    for(g in 1:length(a)) {
print(paste("   ", g))
      sample.subset = which(gen == names(a)[g])
      v_path = matrix(0, length(sample.subset), length(cur.snps))
      v_prob = matrix(0, length(sample.subset), length(cur.snps))

      res = .C(name = viterbi_from_r,
               as.integer(c(length(states), length(sample.subset),
                          length(cur.snps))),
               a = as.double(a[[g]]),
               b = as.double(b[,sample.subset,]),
               v_path = as.integer(v_path),
               v_prob = as.double(v_prob))

    retval[[chr]] = matrix(res$v_path, length(sample.subset), 
                    length(cur.snps), dimnames = 
                    list(rownames(rho)[sample.subset], colnames(rho)[cur.snps]))
    } # for(g)

#  } # for(chr)

  # Next, run the X chromosome for males and females separately.
  

  return(retval)    
} # call.genotypes()
