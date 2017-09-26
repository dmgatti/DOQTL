################################################################################
# Generic skeleton for the alelle call based HMM.
################################################################################
hmm.allele = function(data, founders, sex, snps, chr, trans.prob.fxn) {

  maxIter = 100
  epsilon = 5 # NOTE: We change this in the first iteration below.
  p = 1
  lastLogLik = -Inf
  logLik = -.Machine$double.xmax

  # Convert genotypes to numbers.
  tmp = convert.allele.calls(data$geno, founders$geno)
  data$geno     = tmp[[1]]
  founders$geno = tmp[[2]]
  rm(tmp)

  # Initialize the HMM.
  init.hmm = initialize.hmm(snps = snps[,1], samples = rownames(data$geno),
             states = founders$states)

  # Get emission probabilities.
  b = emission.probs.allele(founders = founders, chr = chr, snps = snps,
      sex = sex)

  # Save the initial emission probabilities as pseudocounts.
  pseudocounts = b
  a = NULL
  if(attr(data, "sampletype") %in% c("DO", "DOF1", "HS", "HSrat")) {

    a = trans.prob.fxn(states = founders$states, snps = snps, chr = chr, 
        sex = sex, gen = data$gen)

  } else {

    a = list(trans.prob.fxn(states = founders$states, snps = snps, chr = chr,
        sex = "F"))

  } # else

  while(p <= maxIter & logLik - lastLogLik > epsilon) {

    print(date())
    print(p)

    # Optional plotting code to watch the HMM progress.
    if(TRUE) {
      prsmth.plot(1, founders$states, init.hmm$prsmth)
    } # if(plot)

    # Initialize the log-likelihood to a large negative number.
    lastLogLik = logLik

    # Filter and smooth each generation separately because they each have a
    # different transition probability matrix.
    print("Running HMM ...")
    lltmp = matrix(-.Machine$double.xmax, 1, 1)

    for(i in 1:length(a)) {

      # Filter
      gen = 1:nrow(data$geno)

      if(any(names(data) == "gen") & attr(data, "sampletype") != "CC") {

        print(paste("gen", names(a)[i]))
        gen = which(data$gen == names(a)[i])

      } # if(any(names(data) == "gen"))

      res = .C(C_filter_smooth_allele,
               dims = as.integer(c(dim(init.hmm$prsmth[,gen,,drop = FALSE]), nrow(b))),
               geno = as.integer(data$geno[gen,]), 
               a = as.double(a[[i]]),
               b = as.double(b),
               prsmth = as.double(init.hmm$prsmth[,gen,,drop = FALSE]),
               init = as.double(init.hmm$init),
               loglik = as.double(lltmp))

       init.hmm$prsmth[,gen,] = res$prsmth
       lltmp = addLog(lltmp, res$loglik)

    } # for(i)

    logLik = lltmp

    print(paste("LogLik =", logLik))

    if(p == 1) {
      epsilon = abs(logLik * 0.001)
    } # if(p == 1)

    # Update the parameters and state means and variances.
    print("Updating Parameters ...")
    b = parameter.update.alleles(geno = data$geno, b = b, 
        pseudocounts = pseudocounts, prsmth = init.hmm$prsmth)
    p = p + 1

  } # while(p < maxIter & logLik - lastLogLik > epsilon

  print(date())

  # If we reached the end of the interations.
  if(p >= maxIter) {
    print("Maximum iterations reached")
  } # if(p >= maxIter)
    
  return(list(b = b, prsmth = init.hmm$prsmth))

} # hmm.allele
