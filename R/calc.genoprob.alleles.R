################################################################################
# Main function for running the haplotype HMM for DO mice.
# You must supply either a set of founders with your data or a set of 
# founder means and variances.  This means that either is.founder.F1 must have
# true values for at least one sample from each founder or init.means and 
# init.covars must be populated with state mean and covariance estimates for all
# SNPs.
# Arguments: data: list, containin
#            chr: character vector, with chromosomes to run.  Must match
#                 the chromosome IDs in the snps table.  NULL means
#                 run all.
#            founders: list, containing genotype calls for the founders
#                      and F1s.
#            snps: data.frame, with 4 columns containing SNPs to use. 
#                  SNP IDs, chr, Mb position, cM position in columns 1:4.
#            output.dir: character, with the file directory for the final
#                        theta & rho mean and variance files.
#            trans.prob: list, containing transition probabilities between
#                        each SNP. Only required for custom sample types.
#            plot: boolean that is true if the user would like to plot a sample
#                  chromosome as the model progresses.
calc.genoprob.alleles = function(data, chr, founders, snps, output.dir = ".",
                        trans.prob.fxn = do.trans.probs, plot = F) {

  # No NAs in geno.
  if(any(is.na(data$geno))) {
    data$geno[is.na(data$geno)] = "N"
    message("Changing NA values in genotypes to 'N'.")
  } # if(any(is.na(data$geno))

  # Select the chromosomes to run.  If chr.to.run = "all", then run all chrs.
  unique.chrs = unique(snps[,2])
  if(chr[1] == "all") {
    chr = unique.chrs
  } # else

  # Loop through each chromosome.
  for(curr.chr in chr) {

    print(paste("CHR", curr.chr))

    # SNPs on the current chromosome.
    cur.snps = which(snps[,2] == curr.chr)

    # If this is the X chromosome, split the samples by sex, using 36 states
    # for the females and 8 states for the males.
    if(curr.chr == "X") {

      # Run the females first.
      females = which(data$sex == "F")
      female.prsmth = NULL
      female.b = NULL

      # Only run if there are samples that are female.
      if(length(females) > 0) {
        print("Females")

        cur.data = list(geno = data$geno[females,cur.snps], 
                   sex = data$sex[females], gen = data$gen[females])
        attributes(cur.data) = attributes(data)
        founder.subset = which(founders$sex == "F")
        cur.founders = list(geno = founders$geno[founder.subset,cur.snps],
                       sex = founders$sex[founder.subset],
                       code = founders$code[founder.subset],
                       states = founders$states$X$F)

        tmp = hmm.allele(data = cur.data, founders = cur.founders,
              sex = "F", snps = snps[cur.snps,], chr = curr.chr, 
              trans.prob.fxn = trans.prob.fxn)
        female.prsmth = tmp$prsmth
        female.b = tmp$b
        rm(tmp)

      } # if(length(females) > 0)

      # Run the males, which only have founder states, not F1s because the 
      # males are hemizygous.
      males = which(data$sex == "M")
      male.prsmth = NULL
      male.b = NULL

      # Only run if there are samples that are male.  If only founders are
      # male, there's no point in running this.
      if(length(males) > 0) {
        print("Males")

        cur.data = list(geno = data$geno[males,cur.snps], 
                   sex = data$sex[males], gen = data$gen[males])
        attributes(cur.data) = attributes(data)
        founder.subset = which(founders$sex == "M")
        cur.founders = list(geno = founders$geno[founder.subset,cur.snps],
                       sex = founders$sex[founder.subset],
                       code = founders$code[founder.subset],
                       states = founders$states$X$M)
        cur.founders = keep.homozygotes(cur.founders)

        tmp = hmm.allele(data = cur.data, founders = cur.founders,
              sex = "M", snps = snps[cur.snps,], chr = curr.chr, 
              trans.prob.fxn = trans.prob.fxn)
        male.prsmth = tmp$prsmth
        male.b = tmp$b
        rm(tmp)

      } # if(length(males) > 0)

      # Combine the male and female prsmth data.
      if(length(females) > 0) {
        if(length(males) > 0) {
          prsmth = array(-.Machine$double.xmax, c(length(founders$states$X$F), 
                   nrow(data$geno), length(cur.snps)), dimnames = 
                   list(founders$states$X$F, rownames(data$geno), 
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
      
      # Write out the smoothed probabilities and the emission probabilities.
      write.results(prsmth = prsmth, b = female.b, output.dir = output.dir, chr = curr.chr,
                    all.chr = chr, sex = "F")
      write.results(b = male.b, output.dir = output.dir, chr = curr.chr,
                    all.chr = chr, sex = "M")
					
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
      x.mchr = x[, rng]
      y.mchr = y[, rng]

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
      cur.data = list(geno = data$geno[,cur.snps], sex = data$sex, gen = data$gen)
      cur.data$geno[cur.data$geno == "-"] = "N"
      attributes(cur.data) = attributes(data)
      cur.founders = list(geno = founders$geno[,cur.snps], sex = founders$sex,
                     code = founders$code, states = founders$states$auto)

      tmp = hmm.allele(data = cur.data, founders = cur.founders,
            sex = "F", snps = snps[cur.snps,], chr = curr.chr,
            trans.prob.fxn = trans.prob.fxn)
      prsmth = tmp$prob.mats$prsmth
      b = tmp$b

      # Write out the smoothed probabilities and the founder state meand and 
      # variances.
      write.results(prsmth = tmp$prsmth, b = tmp$b, output.dir = output.dir, 
	                chr = curr.chr, all.chr = chr)
      rm(tmp)
    } # else
    gc()
	
  } # for(chr)

} # calc.genoprob.alleles()

