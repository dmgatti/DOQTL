################################################################################
# Main function for running the haplotype HMM for DO mice.
# You must supply either a set of founders with your data or a set of 
# founder means and variances.  This means that either is.founder.F1 must have
# true values for at least one sample from each founder or init.means and 
# init.covars must be populated with state mean and covariance estimates for all
# SNPs.
# Arguments: data: list containing x & y intensities or genotypes.
#            chr: character vector, with chromosomes to run.  Must match
#                 the chromosome IDs in the snps table.  NULL means
#                 run all.
#            snps: data.frame with SNP IDs, chr, Mb and cM positions. 
#            output.dir: character, with the file directory for the final
#                        theta & rho mean and variance files.
#            trans.prob.fxn: the transition probability function for these 
#                            samples.
#            plot: boolean that is true if the user would like to plot a sample
#                  chromosome as the model progresses.
calc.genoprob.intensity = function(data, chr, founders, snps, output.dir = ".",
                          trans.prob.fxn, plot = FALSE) {

  # Select the chromosomes to run.  If chr = "all", then run all chrs.
  if(chr[1] == "all") {
    chr = names(snps)
  } # else

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
    if(attr(data, "array") != "megamuga") {
      newxy = quantilenorm(x1 = data$x, y1 = data$y, x2 = founders$x,
              y2 = founders$y)
      founders$x = newxy[[1]]
      founders$y = newxy[[2]]
      rm(newxy)
    } # if(attr(data, "array") != "megamuga")
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
            trans.prob.fxn = trans.prob.fxn)
      founders$states = old.states
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
