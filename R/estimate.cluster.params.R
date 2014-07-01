################################################################################
# Arguments: founders: list with theta, rho, sex and codes for the founders.
#                      Should also include the states for the current chr.
#            chr: character containing the current chromosome.
estimate.cluster.params = function(founders, data, chr) {

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
  theta.means = aggregate(founders$theta, list(founders$code), mean, na.rm = TRUE)
  rho.means   = aggregate(founders$rho,   list(founders$code), mean, na.rm = TRUE)
  founder.means = array(0, c(nrow(theta.means), ncol(founders$rho), 2),
                  dimnames = list(theta.means[,1], colnames(founders$rho),
                  c("theta", "rho")))
  founder.means[,,1] = as.matrix(theta.means[,-1])
  founder.means[,,2] = as.matrix(rho.means[,-1])

  # Loop through each SNP and cluster the samples.
  for(s in 1:ncol(data$rho)) {

    if(s %% 10 == 0) {
      print(paste("SNP", s, "of", ncol(founders$rho)))
    } # if(s %% 10 == 0)

    # Cluster the data.
    rt = cbind(c(founders$theta[,s], data$theta[,s]), c(founders$rho[,s], data$rho[,s]))
    rt = rt[(rowSums(is.na(rt)) + rowSums(is.nan(rt)) + rowSums(is.infinite(rt))) == 0,]
    old.warn = options("warn")$warn
    options(warn = -1)
    mc = Mclust(rt, modelNames = "VVI", prior = priorControl())

    # If we have clusters with less than min.prop of the data, then recluster.
    prop.lt.min.prop = mc$parameters$pro < min.prop

    if(any(prop.lt.min.prop)) {

      G = mc$G - sum(prop.lt.min.prop)
      # If mclust can't calculate a model with this number of clusters, then it
      # will have NA in the BIC value for that number of clusters.  We will
      # cluster using the next lowest number of clusters that mclust can
      # calculate.

      if(is.na(mc$BIC[G])) {
        G = max(which(!is.na(mc$BIC[1:G])))
      } # if(is.na(mc$BIC[G]))
      mc = Mclust(rt, G = 1:G, modelmodelNames = "VVI")

    } # if(any(prop.lt.min.prop))
    options("warn" = old.warn)

    # The temporary classification vector for each state. We will populate
    # this with cluster numbers and then copy the cluster means and variances
    # over.
    cls = rep(0, length(founders$states))
    names(cls) = founders$states

    # If we have one cluster, then every sample is in it.
    if(mc$G == 1) {
      cls = rep(1, length(founders$states))
    } else {

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

  } # for(s)

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
