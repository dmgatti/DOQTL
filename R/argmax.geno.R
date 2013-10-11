################################################################################
# Calculate the Viterbi path through the genotype probabilities for one
# chromosome.
# Daniel Gatti
# Dan.Gatti@jax.org
# Oct. 1, 2012
################################################################################
# Arguments: X.file: character, the path to the X intensities file.
#            Y.file: character, the path to the Y intensities file.
#            sex.gen.file: character, the sex and DO outbreeding generation file.
#            mean.files: character vector, the path to the cluster mean files.
#            covar.files: character vector, the path to the cluster covariance
#                         files.
#            chr: character vector of chromosomes to genotype.
#            snps: data.frame with four columns: SNP ID, Chr, Mb position and
#                  cM position. This can be obtained from data(snps).
#            direction: character, on of "forward" or "bakcward"
# Returns: list of matrices with the Viterbi paths for each Chr.
argmax.geno = function(X.file, Y.file, sex.gen.file, mean.files, covar.files,
              chr, snps, direction = c("forward", "backward")) {

  direction = match.arg(direction)
  data(states)

  # Read in the X & Y files and convert them to rho/theta coordinates.
  x = as.matrix(read.delim(X.file))
  y = as.matrix(read.delim(Y.file))
  x = x[,colnames(x) %in% snps[,1]]
  y = y[,colnames(y) %in% snps[,1]]
  theta = (2.0 / pi) * atan2(y, x)
  rho = x + y
  rm(x, y)

  # Read in the sex and DO outbreeding generation.
  sex.gen = read.delim(sex.gen.file)

  retval = as.list(rep(NA, length(chr)))
  names(retval) = chr
  for(i in 1:length(chr)) {

    print(paste("Processing Chr", chr[i], "..."))
    cur.snps = which(snps[,2] == chr[i])

    # From here on, we have to split between males and females.
    # Males and females can run together on the autosomes, but they have
    # different transition probabilities and genotypes on the X Chr.
    # Autosomes.
    if(chr[i] %in% 1:19) {

      # Read in the cluster means and covariances.
      load(mean.files[grep(paste("chr", chr[i],  "\\.", sep = ""), mean.files)])
      load(covar.files[grep(paste("chr", chr[i], "\\.", sep = ""), covar.files)])

      # Get emission probabilities for this chromosome.
      b = get.log.probabilities(theta = theta[,cur.snps], rho = rho[,cur.snps],
        theta.rho.means = theta.rho.means, theta.rho.covars = theta.rho.covars,
        states = dimnames(theta.rho.means)[[1]])

      # Transition probabilities.
      a = create.log.transition.matrices(dimnames(theta.rho.means)[[1]],
          snps[cur.snps,], do.gen = sort(unique(sex.gen[,2])),
          chr = as.character(chr[i]), sex = "F")

      # Initialize probability matrices. 1/64 for homozygotes, 1/32 for hets.
      init = rep(1.0 / 32.0, dim(b)[1])
      init[c(1,9,16,22,27,31,34,36)] = 1.0 / 64.0
      names(init) = dimnames(b)[[1]]
      init = log(init)

      # Get paths.
      if(direction == "forward") {
        tmp = argmax.helper(a, b, init, sex.gen[,2])
        retval[[i]] = matrix(states[tmp], nrow(tmp), ncol(tmp), dimnames =
                      dimnames(tmp))
      } else {

        tmp = argmax.helper.backward(a, b, init, sex.gen[,2])
        retval[[i]] = matrix(states[tmp], nrow(tmp), ncol(tmp), dimnames =
                      dimnames(tmp))
      } # else

    } else {

      # Female X Chr.
      female.paths = NULL
      female.samples = which(sex.gen[,1] == "F")
      if(length(female.samples) > 0) {

        # Read in the cluster means and covariances.
        load(mean.files[grep(paste("chr", chr[i],  "\\.female", sep = ""),
             mean.files)])
        load(covar.files[grep(paste("chr", chr[i], "\\.female", sep = ""),
             covar.files)])

        # Get emission probabilities for this chromosome.
        b = get.log.probabilities(theta = theta[,cur.snps], rho = rho[,cur.snps],
          theta.rho.means = theta.rho.means, theta.rho.covars = theta.rho.covars,
          states = dimnames(theta.rho.means)[[1]])

        # Transition probabilities.
        a = create.log.transition.matrices(dimnames(theta.rho.means)[[1]],
            snps[cur.snps,], do.gen = sort(unique(sex.gen[,2])),
            chr = as.character(chr[i]), sex = "F")

        # Initialize probability matrices. 1/64 for homozygotes, 1/32 for hets.
        init = rep(1.0 / 32.0, dim(b)[1])
        init[c(1,9,16,22,27,31,34,36)] = 1.0 / 64.0
        names(init) = dimnames(b)[[1]]
        init = log(init)

        # Get paths.
        if(direction == "forward") {
          female.paths = argmax.helper(a, b[,female.samples,], init, 
                         sex.gen[female.samples,2])
        } else {
          female.paths = argmax.helper.backward(a, b[,female.samples,],
                         init, sex.gen[female.samples,2])
        } # else

      } # if(length(female.samples) > 0)

      # Male X Chr.
      male.paths = NULL
      male.samples = which(sex.gen[,1] == "M")
      if(length(male.samples) > 0) {

        # Read in the cluster means and covariances.
        load(mean.files[grep(paste("chr", chr[i],  "\\.male", sep = ""),
             mean.files)])
        load(covar.files[grep(paste("chr", chr[i], "\\.male", sep = ""),
             covar.files)])

        # Get emission probabilities for this chromosome.
        b = get.log.probabilities(theta = theta[,cur.snps], rho = rho[,cur.snps],
          theta.rho.means = theta.rho.means, theta.rho.covars = theta.rho.covars,
          states = dimnames(theta.rho.means)[[1]])

        # Transition probabilities.
        a = create.log.transition.matrices(dimnames(theta.rho.means)[[1]],
            snps[cur.snps,], do.gen = sort(unique(sex.gen[,2])),
            chr = as.character(chr[i]), sex = "M")

        # Initialize probability matrices. 1/64 for homozygotes, 1/32 for hets.
        init = rep(1.0 / 8.0, dim(b)[1])
        names(init) = dimnames(b)[[1]]
        init = log(init)

        # Get paths.
        if(direction == "forward") {
          male.paths = argmax.helper(a, b[,male.samples,], init,
                       sex.gen[male.samples,2])
        } else {
          male.paths = argmax.helper.backward(a, b[,male.samples,],
                       init, sex.gen[male.samples,2])
        } # else

      } # if(length(male.samples) > 0)

      # Assemble the male and female paths.
      paths = matrix(0, dim(b)[2], dim(b)[3], dimnames = list(dimnames(b)[[2]],
              dimnames(b)[[3]]))
      paths[female.samples,] = female.paths
      paths[male.samples,]   = male.paths
      retval[[i]] = matrix(states[paths], nrow(paths), ncol(paths), dimnames =
                    dimnames(paths))
    } # else
  } # for(i)

  return(retval)

} # argmax.geno()


# Helper function to perform the Viterbi algorithm on an arbitrary
# number of states.
argmax.helper = function(a, b, init, do.gen) {

  # t1 hold the maximum probabilities of arriving at each state.
  t1 = array(data = 0, dim = dim(b), dimnames = dimnames(b))
  # t2 holds the state from which we arrived at each current state.
  t2 = array(data = 0, dim = dim(b), dimnames = dimnames(b))

  # Start location
  t1[,,1] = init + b[,,1]
  t2[,,1] = matrix(apply(t1[,,1], 2, which.max), nrow(t2), ncol(t2),
            byrow = T)
  
  num.states  = dim(b)[1]
  num.samples = dim(b)[2]
  num.snps    = dim(b)[3]

  # Run each DO generation separately because they have different transition
  # probabilities.
  for(g in 1:length(a)) {
    sample.subset = which(do.gen == names(a)[g])
    for(i in 2:num.snps) {
      # Multiply the previous probability by the transition prob.
      x = apply(t1[,sample.subset,i-1], 2, "+", a[[g]][,,i-1])
      # Make a 3D array of the values.
      x = array(x, c(num.states, num.states, length(sample.subset)))
      # Get the maximum probability and multiply by the emission prob.
      t1[,sample.subset,i] = apply(x, 2:3, max) + b[,sample.subset,i]
      # Get the state from which we calculated the maximum prob.
      t2[,sample.subset,i] = apply(x, 2:3, which.max)
    } # for(i)
  } # for(g)

  # Trace back and find the most probable path for all samples.
  path = matrix(0, num.samples, num.snps, dimnames = list(dimnames(b)[[2]],
                dimnames(b)[[3]]))
  path[,ncol(path)] = t2[,,dim(t2)[3]][matrix(apply(t1[,,dim(t1)[3]], 2, 
                         which.max), 1:dim(t2)[2])]
  sel.mat = cbind(rep(0, ncol(t2)), 1:ncol(t2))
  for(i in (num.snps-1):1) {
    sel.mat[,1] = path[,i+1]
    path[,i] = t2[,,i][sel.mat]
  } # for(i)
  path[,1] = path[,2]

  return(path)

} # argmax.helper()


# Function to perform Viterbi starting at the distal end of the Chr 
# instead of the proximal end changes the genotypes.
argmax.helper.backward = function(a, b, init, do.gen) {

  # t1 hold the maximum probabilities of arriving at each state.
  t1 = array(data = 0, dim = dim(b), dimnames = dimnames(b))
  # t2 holds the state from which we arrived at each current state.
  t2 = array(data = 0, dim = dim(b), dimnames = dimnames(b))

  # Start location
  t1[,,dim(t1)[3]] = init + b[,,dim(b)[3]]
  t2[,,dim(t2)[3]] = matrix(apply(t1[,,dim(t1)[3]], 2, which.max), nrow(t2),
                     ncol(t2), byrow = T)
  
  num.states  = dim(b)[1]
  num.samples = dim(b)[2]
  num.snps    = dim(b)[3]

  # Run each DO generation separately because they have different transition
  # probabilities.
  for(g in 1:length(a)) {
    sample.subset = which(do.gen == names(a)[g])
    for(i in (num.snps-1):1) {
      # Multiply the previous probability by the transition prob.
      x = apply(t1[,sample.subset,i+1], 2, "+", a[[g]][,,i])
      # Make a 3D array of the values.
      x = array(x, c(num.states, num.states, length(sample.subset)))
      # Get the maximum probability and multiply by the emission prob.
      t1[,sample.subset,i] = apply(x, 2:3, max) + b[,sample.subset,i]
      # Get the state from which we calculated the maximum prob.
      t2[,sample.subset,i] = apply(x, 2:3, which.max)
    } # for(i)
  } # for(g)

  # Trace back and find the most probable path for all samples.
  path = matrix(0, num.samples, num.snps, dimnames = list(dimnames(b)[[2]],
                dimnames(b)[[3]]))
  path[,1] = t2[,,1][matrix(apply(t1[,,1], 2, which.max), 1:dim(t2)[2])]
  sel.mat = cbind(rep(0, ncol(t2)), 1:ncol(t2))
  for(i in 2:num.snps) {
    sel.mat[,1] = path[,i-1]
    path[,i] = t2[,,i][sel.mat]
  } # for(i)
  path[,ncol(path)] = path[,ncol(path)-1]

  return(path)

} # argmax.helper.backward()



