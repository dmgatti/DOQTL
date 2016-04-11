################################################################################
# Arguments: founders: list with theta, rho, sex and codes for the founders.
#                      Should also include the states for the current chr.
#            chr: character containing the current chromosome.
estimate.cluster.params = function(founders, data, chr, 
                                   clust = c("mclust", "pamk")) {

  clust = match.arg(clust)
  min.prop = 0.02
 
  # First verify that we have all of the founders and F1s.
  if(!is.na(as.numeric(chr))) {

    if(!all(founders$states %in% founders$code)) {
      stop(paste("estimate.cluster.params: All founders and F1s not found",
           "in founder data. Please verify that you have one example from",
           "each sex for all founders and F1s."))
    } # if(!all(founders$states %in% founders$code))

  } else if(chr == "X") {

    if(all(founders$sex == "F")) {

      if(!all(founders$states %in% founders$code)) {
        stop(paste("estimate.cluster.params: All founders and F1s not found",
             "in founder data. Please verify that you have one example from",
             "each sex for all founders and F1s."))
      } # if(!all(founders$states %in% founders$code))

    } else if(all(founders$sex == "M")) {

      states = paste(founders$states, founders$states, sep = "")

      if(!all(states %in% founders$code)) {

        stop(paste("estimate.cluster.params: All founders and F1s not found",
             "in founder data. Please verify that you have one example from",
             "each sex for all founders and F1s."))

      } # if(!all(founders$states %in% codes))

    } else {

      stop("estimate.cluster.params: Mixed sex on Chr X not allowed.")

    } # else

  } else if(chr == "Y") {

    if(any(founders$sex == "F")) {

      stop("estimate.cluster.params: Females not allowed on Chr Y.")

    } else {

      states = paste(founders$states, founders$states, sep = "")

      if(!all(states %in% founders$code)) {

        stop(paste("estimate.cluster.params: All founders and F1s not found",
             "in founder data. Please verify that you have one example from",
             "each sex for all founders and F1s."))

      } # if(!all(founders$states %in% codes))

    } # else

  } else if(chr == "M") {

    states = paste(founders$states, founders$states, sep = "")

    if(!all(states %in% founders$code)) {

      stop(paste("estimate.cluster.params: All founders and F1s not found",
         "in founder data. Please verify that you have one example from",
         "each sex for all founders and F1s."))

    } # if(!all(founders$states %in% codes))

  } else {

    stop(paste("estimate.cluster.params: Unknown chr", chr))

  } # else

  r.t.means = array(0, c(length(founders$states), 2, ncol(founders$rho)), 
              dimnames = list(founders$states, c("theta", "rho"),
              colnames(founders$rho)))
  r.t.covars = array(0, dim(r.t.means), dimnames = dimnames(r.t.means))

  # Get the parental means for all SNPs.
  print("Calulating founder mean intensities at each marker...")
  theta.means = aggregate(founders$theta, list(founders$code), mean, na.rm = TRUE)
  rho.means   = aggregate(founders$rho,   list(founders$code), mean, na.rm = TRUE)
  founder.means = array(0, c(nrow(theta.means), ncol(founders$rho), 2),
                  dimnames = list(theta.means[,1], colnames(founders$rho),
                  c("theta", "rho")))
  founder.means[,,1] = as.matrix(theta.means[,-1])
  founder.means[,,2] = as.matrix(rho.means[,-1])

  # Loop through each SNP and cluster the samples.
  print("Estimating intensity clusters at each marker...")
  for(s in 1:ncol(data$rho)) {

    if(s %% 10 == 0) {
      print(paste("Marker", s, "of", ncol(founders$rho)))
    } # if(s %% 10 == 0)

    # Cluster the data.
    rt = cbind(c(founders$theta[,s], data$theta[,s]), c(founders$rho[,s], data$rho[,s]))
    rt = rt[(rowSums(is.na(rt)) + rowSums(is.nan(rt)) + rowSums(is.infinite(rt))) == 0,]
    old.warn = options("warn")$warn
    options(warn = -1)

    if(clust == "mclust") {
      
      mc = Mclust(rt, modelNames = "VVI", prior = priorControl())
      options("warn" = old.warn)

      # The temporary classification vector for each state. We will populate
      # this with cluster numbers and then copy the cluster means and variances
      # over.
      cls = rep(1, length(founders$states))
      names(cls) = founders$states

      # If we have one cluster, then every sample is in it.
      if(mc$G > 1) {

        # Assign founders/F1s to the nearest cluster.
        cl.dist = outer(mc$parameters$mean[1,], founder.means[,s,1], "-")^2 + 
                  outer(mc$parameters$mean[2,], founder.means[,s,2], "-")^2

        cls = apply(cl.dist, 2, which.min)

        if(attr(data, "sampletype") == "DOF1") {
          # We need this for the DO/F1 code to work. It reduces the genotypes
          # classes down to the only possible ones.
          cls = cls[names(cls) %in% founders$states]
        } # if(attr(data, "sampletype") == "DOF1")

      } # else

      # Now, with the cls vector populated, fill in the state means.
      r.t.means[,1,s] = mc$parameters$mean[1,cls]
      r.t.means[,2,s] = mc$parameters$mean[2,cls]

      # Save the variances.
      r.t.covars[,1,s] = mc$parameters$variance$sigma[1,1,cls]
      r.t.covars[,2,s] = mc$parameters$variance$sigma[2,2,cls]
    
    } else {

      mc = pamk(data = rt, krange = 2:9, criterion = "multiasw", usepam = FALSE)
      
      # The temporary classification vector for each state. We will populate
      # this with cluster numbers and then copy the cluster means and variances
      # over.
      cls = rep(1, length(founders$states))
      names(cls) = founders$states

      tmp = split(data.frame(rt), mc$pamobject$clustering)
      tmp = lapply(tmp, as.matrix)
      clmeans = sapply(tmp, colMeans, na.rm = TRUE)
      clcov = sapply(tmp, cov, use = "pairwise.complete.obs")
      clcov = array(clcov, c(2, 2, ncol(clcov)))
      
      # If we have one cluster, then every sample is in it.
      if(mc$nc > 1) {
        
        # Assign founders/F1s to the nearest cluster.
        cl.dist = outer(clmeans[1,], founder.means[,s,1], "-")^2 + 
                  outer(clmeans[2,], founder.means[,s,2], "-")^2
        
        cls = apply(cl.dist, 2, which.min)
        
        if(attr(data, "sampletype") == "DOF1") {
          # We need this for the DO/F1 code to work. It reduces the genotypes
          # classes down to the only possible ones.
          cls = cls[names(cls) %in% founders$states]
        } # if(attr(data, "sampletype") == "DOF1")
        
      } # else
      
      # Now, with the cls vector populated, fill in the state means.
      r.t.means[,1,s] = clmeans[1,cls]
      r.t.means[,2,s] = clmeans[2,cls]
      
      
      # Save the variances.
      r.t.covars[,1,s] = clcov[1,1,cls]
      r.t.covars[,2,s] = clcov[2,2,cls]

    } # else

  } # for(s)

  r.t.covars[r.t.covars < 1e-8] = 1e-8
  
  return(list(r.t.means  = r.t.means, r.t.covars = r.t.covars))

} # estimate.cluster.params()


# Helper function to keep only the homozygote founders.
keep.homozygotes = function(founders) {
  code = matrix(unlist(strsplit(founders$code, split = "")), nrow = 2)
  keep = which(code[1,] == code[2,])
  for(i in 1:length(founders)) {
    if(names(founders)[i] != "states") {
      if(is.matrix(founders[[i]])) {
        founders[[i]] = founders[[i]][keep,]
      } else {
        founders[[i]] = founders[[i]][keep]
      } # else
    } # if(names(founders)[i] != "states")
  } # for(i)
  return(founders)
} # keep.homozygotes()



estimate.cluster.params2 = function(founders, data, chr) {

  min.prop = 0.02
 
  # First verify that we have all of the founders and F1s.
  if(!is.na(as.numeric(chr))) {

    if(!all(founders$states %in% founders$code)) {
      stop(paste("estimate.cluster.params: All founders and F1s not found",
           "in founder data. Please verify that you have one example from",
           "each sex for all founders and F1s."))
    } # if(!all(founders$states %in% founders$code))

  } else if(chr == "X") {

    if(all(founders$sex == "F")) {

      if(!all(founders$states %in% founders$code)) {
        stop(paste("estimate.cluster.params: All founders and F1s not found",
             "in founder data. Please verify that you have one example from",
             "each sex for all founders and F1s."))
      } # if(!all(founders$states %in% founders$code))

    } else if(all(founders$sex == "M")) {

      states = paste(founders$states, founders$states, sep = "")

      if(!all(states %in% founders$code)) {

        stop(paste("estimate.cluster.params: All founders and F1s not found",
             "in founder data. Please verify that you have one example from",
             "each sex for all founders and F1s."))

      } # if(!all(founders$states %in% codes))

    } else {

      stop("estimate.cluster.params: Mixed sex on Chr X not allowed.")

    } # else

  } else if(chr == "Y") {

    if(any(founders$sex == "F")) {

      stop("estimate.cluster.params: Females not allowed on Chr Y.")

    } else {

      states = paste(founders$states, founders$states, sep = "")

      if(!all(states %in% founders$code)) {

        stop(paste("estimate.cluster.params: All founders and F1s not found",
             "in founder data. Please verify that you have one example from",
             "each sex for all founders and F1s."))

      } # if(!all(founders$states %in% codes))

    } # else

  } else if(chr == "M") {

    states = paste(founders$states, founders$states, sep = "")

    if(!all(states %in% founders$code)) {

      stop(paste("estimate.cluster.params: All founders and F1s not found",
         "in founder data. Please verify that you have one example from",
         "each sex for all founders and F1s."))

    } # if(!all(founders$states %in% codes))

  } else {

    stop(paste("estimate.cluster.params: Unknown chr", chr))

  } # else

  x.y.means = array(0, c(length(founders$states), 2, ncol(founders$x)), 
              dimnames = list(founders$states, c("x", "y"),
              colnames(founders$x)))
  x.y.covars = array(0, c(nrow(x.y.means), 3, dim(x.y.means)[[3]]),
               dimnames = list(rownames(x.y.means), c("x", "y", "cov"),
               dimnames(x.y.means)[[3]]))

  # Get the parental means for all SNPs.
  print("Calulating founder mean intensities at each marker...")
  x.means = aggregate(founders$x, list(founders$code), mean, na.rm = TRUE)
  y.means = aggregate(founders$y, list(founders$code), mean, na.rm = TRUE)
  founder.means = array(0, c(nrow(x.means), ncol(founders$y), 2),
                  dimnames = list(x.means[,1], colnames(founders$y),
                  c("x", "y")))
  founder.means[,,1] = as.matrix(x.means[,-1])
  founder.means[,,2] = as.matrix(y.means[,-1])

  # Loop through each SNP and cluster the samples.
  print("Estimating intensity clusters at each marker...")
  for(s in 1:ncol(data$y)) {

    if(s %% 10 == 0) {
      print(paste("Marker", s, "of", ncol(founders$y)))
    } # if(s %% 10 == 0)

    # Cluster the data.
    xy = cbind(c(founders$x[,s], data$x[,s]), c(founders$y[,s], data$y[,s]))
    xy = xy[(rowSums(is.na(xy)) + rowSums(is.nan(xy)) + rowSums(is.infinite(xy))) == 0,]
    old.warn = options("warn")$warn
    options(warn = -1)
    pk = pamk(data = xy, krange = 2:9, criterion = "multiasw", usepam = FALSE)

    options("warn" = old.warn)

    # The temporary classification vector for each state. We will populate
    # this with cluster numbers and then copy the cluster means and variances
    # over.
    cls = rep(1, length(founders$states))
    names(cls) = founders$states

    # Calculate the mean and covariance matrix for each cluster.
    tmp = split(data.frame(xy), pk$pamobject$clustering)
    tmp = lapply(tmp, as.matrix)
    clmeans = sapply(tmp, colMeans, na.rm = TRUE)
    clcov = sapply(tmp, cov, use = "pairwise.complete.obs")

    # If we have more than one cluster, then assign founders to clusters.
    if(pk$nc > 1) {

      # Assign founders/F1s to the nearest cluster.
      cl.dist = outer(clmeans[1,], founder.means[,s,1], "-")^2 + 
                outer(clmeans[2,], founder.means[,s,2], "-")^2

      cls = apply(cl.dist, 2, which.min)

      if(attr(data, "sampletype") == "DOF1") {
        # We need this for the DO/F1 code to work. It reduces the genotypes
        # classes down to the only possible ones.
        cls = cls[names(cls) %in% founders$states]
      } # if(attr(data, "sampletype") == "DOF1")

    } # else

    # Now, with the cls vector populated, fill in the state means.
    x.y.means[,1,s] = clmeans[1,cls]
    x.y.means[,2,s] = clmeans[2,cls]

    # Save the covariances.
    x.y.covars[,1,s] = clcov[1,cls] # X variance in row 1
    x.y.covars[,2,s] = clcov[4,cls] # Y variance in row 4
    x.y.covars[,3,s] = clcov[2,cls] # covariance in row 2 (& 3)

  } # for(s)

  return(list(x.y.means  = x.y.means, x.y.covars = x.y.covars))

} # estimate.cluster.params2()
