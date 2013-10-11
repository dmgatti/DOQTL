################################################################################
# Simulate a QTL with varying effect sizes and splits. We simulate a normally
# distributed trait with mean 0, s.d. 1, add a polygenic effect at various 
# loci and then simulate one, QTL with a specific effect.  We simulate a single
# additive effect. We simulate a single additive QTL on the autosomes. 
# We do not simulate multiple QTL, recessive traits or sex linked traints 
# at this time.
# 
# Daniel Gatti
# Dan.Gatti@jax.org
# Jan. 10, 2013
################################################################################
# Arguments: probs: 3D numeric array, with the founder allele probabilities
#                   as produced by condense.founder.probs(). All dimensions
#                   must have names. Dim 1 contains samples, dim 2 contains
#                   the 8 founders and dim 3 contains SNPs.
#            snps: data.frame with SNP names and locations.  Get from data(snps)
#            K: numeric matrix, with kinship coefficients derived from
#               kinship.probs().  Optional.  If missing, it is ignored.
#            sample.size: numeric vector, with the number of samples to select
#                         from the probs array.  The max must be less than
#                         the number of samples in probs (dim 1).
#            effect.size: numeric vector, with effect sizes to simulate for the
#                         main effect.
#            split: numeric, the number of strains on one side of the allele 
#                   effect split.  For example, 1 would mean an allele split of
#                   1/7.  4 would mean an allele split of 4/4.
#            num.poly: numeric, the number of small polygenic effects to add
#                      to the trait.
#            num.sims: numeric, the number of simulations to perform.
simulate.qtl = function(probs, snps, K, sample.size = dim(probs)[[1]],
               effect.size = 1, split = 4, num.poly = 20, num.sims = 1000) {

  if(missing(probs)) {
    stop(paste("probs cannot be left out from the arguments. Please load",
         "in founder probabilities to pass in."))
  } # if(missing(probs))

  if(max(sample.size) > dim(probs)[[1]]) {
    stop(paste("The maximum number of samples to use (", max(sample.size),
         "cannot exceed the number of samples in the probs array (",
         dim(probs)[[1]],")"))
  } # if(max(sample.size) > dim(probs)[[1]])
 
  probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]
  snps  = snps[snps[,1] %in% dimnames(probs)[[3]],]

  # Create a results list.  We will have one list for each sample size
  # and a list within each sample size for each effect size.  Each effect
  # size will be a matrix with the maximum QTL location, LRS and width.
  results = list()

  # Sample size.
  for(s in 1:length(sample.size)) {
    print(paste("Sample size:", sample.size[s]))
    results[[s]] = list()
    names(results)[s] = sample.size[s]
    # Effect size ...
    for(e in 1:length(effect.size)) {
      print(paste("Effect size:", effect.size[e]))
      results[[s]][[e]] = list(qtl = data.frame(chr = rep("", num.sims),
                          pos = rep(0, num.sims), stringsAsFactors = F),
                          sims = data.frame(chr = rep("", num.sims),
                          pos = rep(0, num.sims), lod = rep(0, num.sims),
                          prox = rep(0, num.sims), dist = rep(0, num.sims),
                          stringsAsFactors = F))
      names(results[[s]])[e] = effect.size[e]
      # Simulation
      for(i in 1:num.sims) {
        print(paste("Simulation:", i))

        # Create the phenotype.
        pheno = rnorm(sample.size[s])

        # Subset the founder probabilities.
        keep = sample(1:dim(probs)[[1]], sample.size[s])
        probs.ss = probs[keep,,]
        if(!missing(K)) {
          K.ss = K[keep, keep]
        } # if(!missing)

        # Simulate a bunch of effects, pick their locations and pick
        # their founder splits.
        loc = sample(1:nrow(snps), num.poly)
        eff = rexp(num.poly)
        eff = eff / sum(eff)
        splits = round(runif(num.poly, 0.5, 4.5))

        # Pick one to be the QTL we seek and set it's effect.
        sim.qtl = sample(1:length(eff), 1)
        eff[sim.qtl] = effect.size[e]
        splits[sim.qtl] = split
        results[[s]][[e]]$qtl[i,] = snps[loc[sim.qtl],2:3]

        # Place the effects among the samples.
        for(j in 1:length(loc)) {
          gt = probs.ss[,,loc[j]]
          strains = sample(1:8, splits[j])
          vec = rep(0, 8)
          vec[strains] = eff[j]
          pheno = pheno + vec %*% t(gt)
        } # for(j)

        # QTL mapping.
        if(missing(K)) {
          qtl = qtl.LRS(pheno = matrix(pheno, ncol = 1), prob = probs.ss,
                snps = snps)$lrs[,c(1:3,5)]
          qtl[,4] = qtl[,4] / (2 * log(10))
        } else {
          qtl = qtl.qtlrel(pheno = matrix(pheno, ncol = 1), prob = probs.ss,
                K = K.ss, snps = snps)$lrs[,c(1:3,7)]
        } # else

        # Record the maximum QTL.
        max.qtl = qtl[which.max(qtl[,4]),]
        interval = bayesint(qtl, chr = max.qtl[1,2])
        results[[s]][[e]]$sims[i,] = c(interval[2,2:3], max.qtl[1,4],
                                       interval[c(1,3),3])
      } # for(i)
    } # for(e)
  } # for(s)

  return(results)

} # simulate.qtl()


