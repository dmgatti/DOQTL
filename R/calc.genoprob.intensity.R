################################################################################
# Main function for running the haplotype HMM for DO mice.
# You must supply either a set of founders with your data or a set of 
# founder means and variances.  This means that either is.founder.F1 must have
# true values for at least one sample from each founder or init.means and 
# init.covars must be populated with state mean and covariance estimates for all
# SNPs.
# Arguments: data: list containing x & y intensities or genotypes.
#            chr: character vector, with chromosomes to run.  Must match
#                 the chromosome IDs in the snps table.
#            snps: data.frame with SNP IDs, chr, Mb and cM positions. 
#            output.dir: character, with the file directory for the final
#                        theta & rho mean and variance files.
#            trans.prob.fxn: the transition probability function for these 
#                            samples.
#            plot: boolean that is true if the user would like to plot a sample
#                  chromosome as the model progresses.
calc.genoprob.intensity = function(data, chr, founders, snps, output.dir = ".",
                          trans.prob.fxn, plot = FALSE, clust = c("mclust", "pamk")) {
  clust = match.arg(clust)
  
  # Extract the sample and founder temporary data file names.
  tmpfiles = list(dx = data$x, dy = data$y, fx = founders$x, fy = founders$y)
  # Loop through each chromosome.
  for(curr.chr in chr) {
    print(paste("CHR", curr.chr))
    # Get SNPs for the current chromosome.
    cur.snps = snps[[which(names(snps) == curr.chr)]]
    # Read in the sample and founder data.
    d = NULL # d is the variable that is loaded in the load() statements below.
    chr.pattern = paste("chr", curr.chr, "\\.", sep = "")
    load(tmpfiles$dx[grep(chr.pattern, tmpfiles$dx)])
    data$x = d
    load(tmpfiles$dy[grep(chr.pattern, tmpfiles$dy)])
    data$y = d
    load(tmpfiles$fx[grep(chr.pattern, tmpfiles$fx)])
    founders$x = d
    load(tmpfiles$fy[grep(chr.pattern, tmpfiles$fy)])
    founders$y = d
    rm(d)
    gc()
    # Normalize the samples and founders.
    # MegaMUGA founders seem to be well aligned to the data, so we're
    # not normalizing MegaMUGA data at the moment.
    if(attr(data, "array") == "muga") {
      newxy = quantilenorm(x1 = data$x, y1 = data$y, x2 = founders$x,
              y2 = founders$y)
      founders$x = newxy[[1]]
      founders$y = newxy[[2]]
      rm(newxy)
    } # if(attr(data, "array") != "muga")
    gc()
    # If this is the X chromosome, split the samples by sex, using 36 states
    # for the females and 8 states for the males.
    if(curr.chr == "X") {
      # Run the females first.
      females = which(data$sex == "F")
      female.prsmth = NULL
      female.r.t.means  = NULL
      female.r.t.covars = NULL
      # Only run if there are samples that are female.
      if(length(females) > 0) {
        print("Females")
        data$theta = (2.0 / pi) * atan2(data$y, data$x)
        data$rho = sqrt(data$x^2 + data$y^2)
        cur.data = list(theta = (2.0 / pi) * atan2(data$y[females,], 
                   data$x[females,]), rho = sqrt(data$x[females,]^2 + 
                   data$y[females,]^2), 
                   sex = data$sex[females], gen = data$gen[females])
        attr(cur.data, "sampletype") = attr(data, "sampletype")
        founder.subset = which(founders$sex == "F")
        cur.founders = list(theta = (2.0 / pi) * atan2(founders$y[founder.subset,], 
                       founders$x[founder.subset,]), 
                       rho = sqrt(founders$x[founder.subset,]^2 + 
                       founders$y[founder.subset,]^2), 
                       sex = founders$sex[founder.subset],
                       code = founders$code[founder.subset],
                       states = founders$states$X$F)
        tmp = hmm.intensity(data = cur.data, founders = cur.founders,
              sex = "F", snps = cur.snps, chr = curr.chr, 
              trans.prob.fxn = trans.prob.fxn, clust = clust)
        female.r.t.means  = tmp$params$r.t.means
        female.r.t.covars = tmp$params$r.t.covars
        female.prsmth = tmp$prsmth
        rm(tmp)
      } # if(length(females) > 0)
      males = which(data$sex == "M")
      male.prsmth = NULL
      male.r.t.means  = NULL
      male.r.t.covars = NULL
      # Only run if there are samples that are male.  If only founders are
      # male, there's no point in running this.
      if(length(males) > 0) {
        print("Males")
        cur.data = list(theta = (2.0 / pi) * atan2(data$y[males,], 
                   data$x[males,]), rho = sqrt(data$x[males,]^2 + 
                   data$y[males,]^2), sex = data$sex[males],
                   gen = data$gen[males])
        attr(cur.data, "sampletype") = attr(data, "sampletype")
        founder.subset = which(founders$sex == "M")
        cur.founders = list(theta = (2.0 / pi) * atan2(founders$y[founder.subset,], 
                       founders$x[founder.subset,]), 
                       rho = sqrt(founders$x[founder.subset,]^2 + 
                       founders$y[founder.subset,]^2), 
                       sex = founders$sex[founder.subset],
                       code = founders$code[founder.subset],
                       states = founders$states$X$M)
        cur.founders = keep.homozygotes(cur.founders)
        tmp = hmm.intensity(data = cur.data, founders = cur.founders,
              sex = "M", snps = cur.snps, chr = curr.chr,
              trans.prob.fxn = trans.prob.fxn, clust = clust)
        male.r.t.means  = tmp$params$r.t.means
        male.r.t.covars = tmp$params$r.t.covars
        male.prsmth = tmp$prsmth
        rm(tmp)
      } # if(length(males) > 0)
	  
      # Combine the male and female prsmth data.
      if(length(females) > 0) {
        if(length(males) > 0) {
          prsmth = array(-.Machine$double.xmax, c(length(founders$states$X$F), 
                   nrow(data$x), nrow(cur.snps)), dimnames = 
                   list(founders$states$X$F, rownames(data$x), cur.snps[,1]))
          m = match(dimnames(female.prsmth)[[2]], dimnames(prsmth)[[2]])
          prsmth[,m,] = female.prsmth
          m = match(dimnames(male.prsmth)[[2]], dimnames(prsmth)[[2]])
          m2 = match(paste(dimnames(male.prsmth)[[1]], dimnames(male.prsmth)[[1]], 
               sep = ""), dimnames(prsmth)[[1]])
          prsmth[m2,m,] = male.prsmth
          rm(female.prsmth, male.prsmth)
        } else {
         prsmth = female.prsmth
         rm(female.prsmth)
        } # else
      } else if(length(males) > 0) {
         prsmth = male.prsmth
         rm(male.prsmth)
      } # else if(length(males) > 0)
      
      # Write out the smoothed probabilities and the founder state means and 
      # variances.
      if(!is.null(female.r.t.means)) {
        write.results(prsmth = prsmth, theta.rho.means = female.r.t.means, 
                      theta.rho.covars = female.r.t.covars,
                      output.dir = output.dir, chr = curr.chr, all.chr = chr,
                      sex = "F")
      } # if(!is.null(female.r.t.means))
      if(!is.null(male.r.t.means)) {
        # We don't want to pass prsmth in twice or else it will be written
        # out twice.
        write.results(theta.rho.means = male.r.t.means, theta.rho.covars = 
                      male.r.t.covars, output.dir = output.dir, chr = curr.chr, 
                      all.chr = chr, sex = "M")
      } # if(!is.null(male.r.t.means))
    } else if (curr.chr == "Y") {
      stop("Chr Y not implemented yet.")
### DMG: NOT IMPLEMENTED YET.
      # Keep only the males.
      males = which(data$sex == "M")
      founder.males = which(founders$sex == "M")
      # Get founder centers.
      fx = founders$x[founder.males,]
      fy = founders$y[founder.males,]
      fx = aggregate(fx, list(founders$code[founder.males]), mean)
      fy = aggregate(fy, list(founders$code[founder.males]), mean)
      # Get distances between founders and samples.
      x = t(rbind(founders$x[founder.males,], data$x[males,]))
      y = t(rbind(founders$y[founder.males,], data$y[males,]))
      d = matrix(0, ncol(x), ncol(x), dimnames = list(colnames(x), colnames(x)))
      for(i in 1:ncol(x)) {
        d[i,] = sqrt(colSums((x[,i] - x)^2) + colSums((y[,i] - y)^2))
      } # for(i)
      d = d[!is.na(d[1,]), !is.na(d[,1])]
      # Cluster using CLARA.
#      cl = pamk(x = d, usepam = FALSE)
  
      # Calculate means and covariances for each cluster.
#      mean.covar = as.list(1:length(founders$states$founders))
#      names(mean.covar) = founders$states$founders
#      for(i in 1:length(founders$states$founders)) {
#        rng = which(cl$clustering == i)
#        mean.covar[[i]]$mean = cbind(rowMeans(x[,rng]), rowMeans(y[,rng]))
#        mean.covar[[i]]$covar = cov(t(x[,rng]), t(y[,rng]))
#      } # for(i)
      
      
    } else if (curr.chr == "M") {
      stop("Chr M not implemented yet.")
    } else {
      # Autosomes.
      data$theta = (2.0 / pi) * atan2(data$y, data$x)
      data$rho = sqrt(data$x^2 + data$y^2)
      founders$theta = (2.0 / pi) * atan2(founders$y, founders$x)
      founders$rho = sqrt(founders$x^2 + founders$y^2)
      old.states = founders$states
      founders$states = founders$states$auto
      hmm = hmm.intensity(data = data, founders = founders,
            sex = "F", snps = cur.snps, chr = curr.chr,
            trans.prob.fxn = trans.prob.fxn, clust = clust)
      founders$states = old.states
      # Write out the smoothed probabilities and the founder state means and 
      # variances.
      write.results(prsmth = hmm$prsmth, 
                    theta.rho.means  = hmm$params$r.t.means, 
                    theta.rho.covars = hmm$params$r.t.covars, 
                    output.dir = output.dir, chr = curr.chr, all.chr = chr)
    } # else
    # Clean up memory.
    gc()
    # Make sure the foreach doesn't try to return something.
    NULL
  } # for(chr)
} # calc.genoprob.intensity()
