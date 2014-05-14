################################################################################
# Given a smoothed probability *.Rdata file, read it in and determine the
# founder allele breakpoints by splitting the change location in half.
# Daniel Gatti
# Dan.Gatti@jax.org.
# Dec. 4, 2012
################################################################################
# Arguments: filenames: character vector with *.genotype.probs.Rdata filenames.
#                       By default, the function looks for 
#                       *.genotype.probs.Rdata files in the current directory.
#                       snps: data.frame containing SNO ID, chr, Mb and cM
#                             for each marker in columns 1:4, respectively.
write.founder.genomes = function(filenames = dir(path = ".",
                        pattern = "genotype.probs.Rdata"), snps) {
  if(length(filenames) > 0) {
    load(filenames[1])
    snps = snps[snps[,1] %in% dimnames(prsmth)[[1]],]
    
    chrlen = get.chr.lengths()
    for(i in 1:length(filenames)) {
      print(filenames[i])
      load(filenames[i]) # this loads prsmth.
      prsmth = prsmth[snps[,1],]
      # Get the maximum genotype at each SNP.
      max.geno = colnames(prsmth)[apply(prsmth, 1, which.max)]
      # Split the genotypes up by chormosome.
      max.geno = split(max.geno, snps[,2])
      max.geno = max.geno[c(1, 12:19, 2:11, 20)]
      # For each chromosome, make breakpoint lists for each strand.
      strand1 = rep(NA, 3)
      strand2 = rep(NA, 3)
      for(c in 1:length(max.geno)) {
        # Split up the two letter genotypes codes.
        max.geno[[c]] = strsplit(max.geno[[c]], split = "")
        max.geno[[c]] = matrix(unlist(max.geno[[c]]), ncol = 2, byrow = TRUE)
        max.geno[[c]] = cbind(snps[snps[,2] == names(max.geno)[c], 2:3],
                        max.geno[[c]])
        max.geno[[c]][,3] = as.character(max.geno[[c]][,3])
        max.geno[[c]][,4] = as.character(max.geno[[c]][,4])
        # Try to minimize recombinations.
        for(j in 2:nrow(max.geno[[c]])) {
          if(max.geno[[c]][j-1,3] == max.geno[[c]][j,4] | 
             max.geno[[c]][j,3]   == max.geno[[c]][j-1,4]) {
            swap = max.geno[[c]][j,3]
            max.geno[[c]][j,3] = max.geno[[c]][j,4]
            max.geno[[c]][j,4] = swap
          } # if(max.geno[[c]][j-1,3] == max.geno[[c]][j,4] |  ...
        } # for(j)
        # Get breakpoints on each strand.
        b1 = which(max.geno[[c]][-1,3] != max.geno[[c]][-nrow(max.geno[[c]]),3])
        b2 = which(max.geno[[c]][-1,4] != max.geno[[c]][-nrow(max.geno[[c]]),4])
        # Strand 1
        s1 = matrix(0, length(b1) + 1, 3, dimnames = list(NULL, c("Chr",
                  "End", "1")))
        s1[1,2:3] = c(1e-6, max.geno[[c]][1,3])
        s1[,1] = max.geno[[c]][1,1]
        s1[-1,2] = 0.5 * (max.geno[[c]][b1,2] + max.geno[[c]][b1+1,2])
        s1[,3] = c(max.geno[[c]][b1,3], max.geno[[c]][nrow(max.geno[[c]]),3])
        s1[,2] = format(round(as.numeric(s1[,2]) * 1e6), digits = 6)
        # Strand 2        
        s2 = matrix(0, length(b2) + 1, 3, dimnames = list(NULL, c("Chr",
                  "End", "2")))
        s2[1,2:3] = c(1e-6, max.geno[[c]][1,4])
        s2[,1] = max.geno[[c]][1,1]
        s2[-1,2] = 0.5 * (max.geno[[c]][b2,2] + max.geno[[c]][b2+1,2])
        s2[,3] = c(max.geno[[c]][b2,4], max.geno[[c]][nrow(max.geno[[c]]),4])
        s2[,2] = format(round(as.numeric(s2[,2]) * 1e6), digits = 6)
        strand1 = rbind(strand1, s1)
        strand2 = rbind(strand2, s2)
      } # for(c)
      strand1 = strand1[-1,]
      strand2 = strand2[-1,]
      sampleID = sub("\\.genotype\\.probs\\.Rdata", "", filenames[i])
      write.csv(strand1, paste(sampleID, ".founder.blocks.A.csv", sep = ""),
                row.names = FALSE, quote = FALSE)
      write.csv(strand2, paste(sampleID, ".founder.blocks.B.csv", sep = ""),
                row.names = FALSE, quote = FALSE)
    } # for(i)
  } # if(length(filenames) > 0)
} # write.founder.genomes()
