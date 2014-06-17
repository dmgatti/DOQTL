################################################################################
# Simulate a QTL with varying effect sizes and mafs. We simulate a normally
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
#            maf: numeric, the number of strains on one side of the allele 
#                   effect maf.  For example, 1 would mean an allele maf of
#                   1/7.  4 would mean an allele maf of 4/4.
#            num.poly: numeric, the number of small polygenic effects to add
#                      to the trait.
#            num.sims: numeric, the number of simulations to perform.
qtl.simulate = function(probs, snps, K, sample.size = dim(probs)[[1]],
               effect.size = 1, maf = 4, num.poly = 18, num.sims = 1000) {
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
      qtl = as.list(1:num.sims)
      names(qtl) = 1:num.sims
      pheno = list(1:num.sims)
	  names(pheno) = 1:num.sims
      # Simulation
      for(i in 1:num.sims) {
        if(i %%100 == 0) print(paste("Simulation:", i))
        # Create the QTL matrix.
        qtl[[i]] = data.frame(chr = 1:19, loc = rep(0, 19), effect = rep(0, 19),
                   maf = rep("", 19), var.expl = rep(0, 19), sim = rep(FALSE, 19))
        # Select the sample subset.
        ss = sort(sample(1:dim(probs)[1], s))
        # Simulate noise with mean = 0, and Std. Dev. = 1.0.
        ph = rnorm(n = s, mean = 0, sd = 1.0)
        names(ph) = dimnames(probs)[[1]][ss]
        # Subset the founder probabilities.
        probs.ss = probs[ss,,]
        # Simulate effects, pick their locations and pick
        # their founder mafs.
        loc = 1:19
        for(j in 1:length(loc)) {
          curr.chr = which(snps[,2] == j)
          loc[j]   = sample(curr.chr, 1)
        } # for(j)
        # Simulate exponentially distributed effects.
        eff = rexp(19)
        mafs = round(runif(19, 0.5, 4.5))
        # Pick one to be the QTL we seek and set it's effect.
        sim.qtl = sample(1:length(eff), 1)
        eff[sim.qtl] = e
        mafs[sim.qtl] = maf
        qtl[[i]]$loc = snps[loc,3] * 1e-6
        qtl[[i]]$maf  = mafs
        qtl[[i]]$sim[sim.qtl] = TRUE
        # Place the effects among the samples.
        gen.effect = matrix(0, s, 19)
        gt = matrix(0, s, 19)
        for(j in 1:length(loc)) {
          # Pick the strains with the effect.
          strains = sample(1:8, mafs[j])
          vec = rep(0, 8)
          vec[strains] = 1
          # Round and scale the genotypes to -1, 0 & 1.
          gt[,j] = round((2 * (vec %*% t(probs.ss[,,loc[j]]))[1,]) - 1)
          # Save the effect.
          gen.effect[,j] = gt[,j] * eff[j]
        } # for(j)
        # Get the variances of the polygenic background.
        gen.var = apply(gen.effect[,-sim.qtl], 2, var)
        # Scale the sum of the variances to 1.
        gen.var = gen.var / sum(gen.var)
        # Adjust the genetic effects to sum to 1.
        tmp = gen.effect[,-sim.qtl]
        tmp = (tmp - matrix(colMeans(tmp), nrow(tmp), ncol(tmp),
              byrow = TRUE))^2
        adj = sqrt(colSums(tmp) / (nrow(tmp) - 1) / gen.var)
        gen.effect[,-sim.qtl] = gen.effect[,-sim.qtl] / 
              matrix(adj, nrow(gen.effect), ncol(gen.effect) - 1,
              byrow = TRUE)
        qtl[[i]]$effect = apply(apply(abs(gen.effect), 2, unique), 2, sort)[2,]
#print(paste("error.var. =", var(ph)))
#print(paste("gen.var. =", var(rowSums(gen.effect[,-sim.qtl]))))
        if(is.na(var(rowSums(gen.effect[,-sim.qtl])))) {
          stop("!!!")
        }
        # Add the polygenic effect to the phenotype.
        ph = ph + rowSums(gen.effect[,-sim.qtl])
        # Standardize the phenotype.
        ph = scale(ph)[,1]
        # Add the simulated QTL to the phenotype.
        j = sim.qtl
        strains = sample(1:8, mafs[j])
        vec = rep(0, 8)
        vec[strains] = 1
        gt[,j] = round((2 * (vec %*% t(probs.ss[,,loc[j]]))[1,]) - 1)
        ph = ph + gt[,j] * eff[j]
        pheno[[i]] = ph
        # Fit a model and calculate the % variance explained.
        for(j in 1:length(loc)) {
          qtl[[i]]$var.expl[j] = cor(ph, gt[,j])^2
        } # for(j)
        # Print out the total variance.
#        print(paste("var =", var(ph)))
      } # for(i)
      save(pheno, qtl, file = paste("pheno.qtl.maf", maf, 
           "sample.size", s ,"effect.size", e, "Rdata", sep = "."))
    } # for(e)
  } # for(s)
  return(results)
} # qtl.simulate()
