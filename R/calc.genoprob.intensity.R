################################################################################
# Main function for running the haplotype HMM for DO mice.
# You must supply either a set of founders with your data or a set of 
# founder means and variances.  This means that either is.founder.F1 must have
# true values for at least one sample from each founder or init.means and 
# init.covars must be populated with state mean and covariance estimates for all
# SNPs.
# Arguments: x: matrix, num.samples x num.snps, with x intensities for all
#               samples. Cannot be passed in if geno is used.
#            y: matrix, num.samples x num.snps, with y intensities for all
#               samples. Cannot be passed in if geno is used.
#            geno: matrix, num.samples x num.snps, with genotypes for all
#               samples. Cannot be passed in if x & y are used.
#            sex.gen: data.frame, with two named columns.  'sex', which contains
#                     either M or F for all samples, and 'gen' which contains the
#                     generation of DO outbreeding. Rownames must contain sample IDs.
#            chr: character vector, with chromosomes to run.  Must match
#                 the chromosome IDs in the snps table.  NULL means
#                 run all.
#            output.dir: character, with the file directory for the final
#                        theta & rho mean and variance files.
#            plot: boolean that is true if the user would like to plot a sample
#                  chromosome as the model progresses.
calc.genoprob.intensity = function(data, chr, founders, snps, output.dir = ".",
                          trans.prob.fxn, plot = F) {

  # Select the chromosomes to run.  If chr = "all, then run all chrs.
  unique.chrs = unique(snps[,2])
  if(chr[1] == "all") {
    chr = unique.chrs[unique.chrs %in% snps$Chr]
  } # else

  # Normalize the intensity values between the founders and data.
  batch = rep(1:2, c(nrow(founders$x), nrow(data$x)))
  xy = normalize.batches(x = rbind(founders$x, data$x),
                         y = rbind(founders$y, data$y),
                         sex = c(founders$sex, data$sex),
                         batch = batch, snps = snps)
  founders$x = xy$x[batch == 1,]
  founders$y = xy$y[batch == 1,]
  data$x = xy$x[batch == 2,]
  data$y = xy$y[batch == 2,]

  # Loop through each chromosome.
  for(curr.chr in chr) {
    print(paste("CHR", curr.chr))

    # theta & rho values for the current chromosome.
    cur.snps = which(snps[,2] == curr.chr)

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

        cur.data = list(theta = (2.0 / pi) * atan2(data$y[females,cur.snps], 
                   data$x[females,cur.snps]), rho = sqrt(data$x[females,cur.snps]^2 + 
                   data$y[females,cur.snps]^2), 
                   sex = data$sex[females], gen = data$gen[females])
        attributes(cur.data) = attributes(data)
        names(cur.data)[1:2] = c("theta", "rho")
        founder.subset = which(founders$sex == "F")
        cur.founders = list(theta = (2.0 / pi) * atan2(founders$y[founder.subset,cur.snps], 
                       founders$x[founder.subset,cur.snps]), 
                       rho = sqrt(founders$x[founder.subset,cur.snps]^2 + 
                       founders$y[founder.subset,cur.snps]^2), 
                       sex = founders$sex[founder.subset],
                       code = founders$code[founder.subset],
                       states = founders$states$X$F)

        tmp = hmm.intensity(data = cur.data, founders = cur.founders,
              sex = "F", snps = snps[cur.snps,], chr = curr.chr, 
              trans.prob.fxn = trans.prob.fxn)
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

        cur.data = list(theta = (2.0 / pi) * atan2(data$y[males,cur.snps], 
                   data$x[males,cur.snps]), rho = sqrt(data$x[males,cur.snps]^2 + 
                   data$y[males,cur.snps]^2), sex = data$sex[males],
                   gen = data$gen[males])
        attributes(cur.data) = attributes(data)
        names(cur.data)[1:2] = c("theta", "rho")
        founder.subset = which(founders$sex == "M")
        cur.founders = list(theta = (2.0 / pi) * atan2(founders$y[founder.subset,cur.snps], 
                       founders$x[founder.subset,cur.snps]), 
                       rho = sqrt(founders$x[founder.subset,cur.snps]^2 + 
                       founders$y[founder.subset,cur.snps]^2), 
                       sex = founders$sex[founder.subset],
                       code = founders$code[founder.subset],
                       states = founders$states$X$M)
        cur.founders = keep.homozygotes(cur.founders)

        tmp = hmm.intensity(data = cur.data, founders = cur.founders,
              sex = "M", snps = snps[cur.snps,], chr = curr.chr, 
              trans.prob.fxn = trans.prob.fxn)
        male.r.t.means  = tmp$params$r.t.means
        male.r.t.covars = tmp$params$r.t.covars
        male.prsmth = tmp$prsmth
        rm(tmp)

      } # if(length(males) > 0)

      # Combine the male and female prsmth data.
      if(length(females) > 0) {
        if(length(males) > 0) {
          prsmth = array(-.Machine$double.xmax, c(length(founders$states$X$F), 
                   nrow(data$x), length(cur.snps)), dimnames = 
                   list(founders$states$X$F, rownames(data$x), 
                   snps[cur.snps,1]))
          m = match(dimnames(female.prsmth)[[2]], dimnames(prsmth)[[2]])
          prsmth[,m,] = female.prsmth
          m = match(dimnames(male.prsmth)[[2]], dimnames(prsmth)[[2]])
          m2 = match(paste(dimnames(male.prsmth)[[1]], dimnames(male.prsmth)[[1]], sep = ""),
               dimnames(prsmth)[[1]])
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

      # Get the male Y chr snps.
      rng  = which(snps[,2] == "Y")
      male = which(sex == "M")
      x.ychr = x[male, rng]
      y.ychr = y[male, rng]

      # Get the founder centers.
      founder.rows = which(rownames(x.ychr) %in% founder.names)
      sample.rows = setdiff(1:nrow(x.ychr), founder.rows)
      xmean = apply(x.ychr[founder.rows,], 2, tapply,
              rownames(x.ychr)[founder.rows], mean)
      ymean = apply(y.ychr[founder.rows,], 2, tapply,
              rownames(y.ychr)[founder.rows], mean)

      # Get the distance of each sample from the founder means.
      dist = matrix(0, length(sample.rows), nrow(xmean), dimnames = list(
             rownames(x.ychr)[sample.rows], rownames(xmean)))
      for(s in 1:ncol(x.ychr)) {
        dist = dist + outer(x.ychr[sample.rows, s], xmean[,s], "-")^2 +
                      outer(y.ychr[sample.rows, s], ymean[,s], "-")^2
      } # for(s)

      min.dist = apply(dist, 1, which.min)

    } else if (curr.chr == "M") {

      # Get the male Y chr snps.
      rng  = which(snps[,2] == "M")
      x.mchr = data$x[, rng]
      y.mchr = data$y[, rng]

      # Get the founder centers.
      founder.rows = which(rownames(x.mchr) %in% founder.names)
      sample.rows  = setdiff(1:nrow(x.mchr), founder.rows)
      xmean = apply(x.mchr[founder.rows,], 2, tapply,
              rownames(x.mchr)[founder.rows], mean)
      ymean = apply(y.mchr[founder.rows,], 2, tapply,
              rownames(y.mchr)[founder.rows], mean)

      # Get the distance of each sample from the founder means.
      dist = matrix(0, length(sample.rows), nrow(xmean), dimnames = list(
             rownames(x.mchr)[sample.rows], rownames(xmean)))
      for(s in 1:ncol(x.mchr)) {
        dist = dist + outer(x.mchr[sample.rows, s], xmean[,s], "-")^2 +
                      outer(y.mchr[sample.rows, s], ymean[,s], "-")^2
      } # for(s)

      min.dist = apply(dist, 1, which.min)

    } else {

      # Autosomes.
      cur.data = list(theta = (2.0 / pi) * atan2(data$y[,cur.snps], 
                 data$x[,cur.snps]), rho = sqrt(data$x[,cur.snps]^2 + 
                 data$y[,cur.snps]^2), sex = data$sex, gen = data$gen)
      attributes(cur.data) = attributes(data)
      names(cur.data)[1:2] = c("theta", "rho")
      cur.founders = list(theta = (2.0 / pi) * atan2(founders$y[,cur.snps], 
                     founders$x[,cur.snps]), rho = sqrt(founders$x[,cur.snps]^2 + 
                     founders$y[,cur.snps]^2), sex = founders$sex,
                     code = founders$code, states = founders$states$auto)

      hmm = hmm.intensity(data = cur.data, founders = cur.founders,
            sex = "F", snps = snps[cur.snps,], chr = curr.chr,
            trans.prob.fxn = trans.prob.fxn)

      # Write out the smoothed probabilities and the founder state means and 
      # variances.
      write.results(prsmth = hmm$prsmth, 
                    theta.rho.means  = hmm$params$r.t.means, 
                    theta.rho.covars = hmm$params$r.t.covars, 
                    output.dir = output.dir, chr = curr.chr, all.chr = chr)
    } # else

	gc()
  } # for(chr)

} # calc.genoprob.intensity()

