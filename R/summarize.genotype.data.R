################################################################################
# Read through a directory of genotype probabilities from the HMM
# (*.genotype.probs.Rdata) and summarize the data by sample and SNP.
# Note: the *.genotype.probs.Rdata files must have already been created.
# The output is a file with all genotype transitions for each sample.
# Daniel Gatti
# Dan.Gatti@jax.org
# April 16, 2012
################################################################################
summarize.genotype.transitions = function(path = ".", snps) {

  path = add.slash(path)
  files = dir(path = path, pattern = "genotype.probs.Rdata", full.names = T)

  if(length(files) == 0) {
    stop(paste("There are no *.genotype.probs.Rdata files in this directory."))
  } # if(length(files) == 0)

  # Load in one file to get the SNPs.
  load(files[1])
  snps = snps[snps[,1] %in% rownames(prsmth),]
  snp.list = split(snps, snps[,2])
  snp.list = snp.list[order(as.numeric(names(snp.list)))]

  # Create a set of 2,000,000 row vectors.  I tried to use a data.frame
  # at this stage and it was very slow.
  nrow = 2e6
  sample = rep(NA, nrow)
  prox.geno = rep(NA, nrow)
  dist.geno = rep(NA, nrow)
  prox.snp = rep(NA, nrow)
  dist.snp = rep(NA, nrow)
  chr = rep(NA, nrow)
  prox.loc = rep(NA, nrow)
  dist.loc = rep(NA, nrow)

  index = 1
  for(f in files) {
    print(f)
    # This loads 'prsmth' that contains genotype probabilities on a log scale.
    load(f)
    prsmth = exp(prsmth)
    prsmth = prsmth[rownames(prsmth) %in% snps[,1],]
    prsmth = split(data.frame(prsmth), snps[,2])
    prsmth = prsmth[order(as.numeric(names(prsmth)))]
    prsmth = lapply(prsmth, as.matrix)

    for(c in 1:length(prsmth)) {
      mx.state = colnames(prsmth[[c]])[apply(prsmth[[c]], 1, which.max)]
      chg = which(mx.state[-1] != mx.state[-length(mx.state)])
      if(length(chg) > 0) {
        rng = index:(index + length(chg) - 1)
        sample[rng] = gsub(paste(path, "|\\.genotype\\.probs\\.Rdata",
                                   sep = ""), "", f)
        prox.geno[rng] = mx.state[chg]
        dist.geno[rng] = mx.state[chg+1]
        chr[rng] = names(prsmth)[c]
        prox.snp[rng] = snp.list[[c]][chg,1]
        dist.snp[rng] = snp.list[[c]][chg+1,1]
        prox.loc[rng] = snp.list[[c]][chg,3]
        dist.loc[rng] = snp.list[[c]][chg+1,3]
        index = index + length(chg)
      } # if(length(chg) > 0)
    } # for(c)
  } # for(f)

  # Remove the unused rows.
  remove = which(is.na(sample))
  sample = sample[-remove]
  prox.geno = prox.geno[-remove]
  dist.geno = dist.geno[-remove]
  prox.snp = prox.snp[-remove]
  dist.snp = dist.snp[-remove]
  chr  = chr[-remove]
  prox.loc = prox.loc[-remove]
  dist.loc = dist.loc[-remove]

  # Create a data.frame and return.
  return(data.frame(sample, prox.geno, dist.geno, prox.snp, dist.snp, 
         chr, prox.loc, dist.loc))
} # summarize.genotype.transitions()


################################################################################
# Go through the genotype probability files and summarize by SNP.
# Arguments: path: character, full path to the directory with the
#                  *.genotype.probs.Rdata files.
#            snps: data.frame with 4 columns containing the SNP ID, Chr, Mb &
#                  cM position of each SNP on the genotyping array.
# Returns: data.frame with one row per SNP and one column for each genotype
#          and founder.
summarize.by.snps = function(path = ".", snps) {

  path = add.slash(path)
  files = dir(path = path, pattern = "genotype.probs.Rdata", full.names = T)

  if(length(files) == 0) {
    stop(paste("There are no *.genotype.probs.Rdata files in this directory."))
  } # if(length(files) == 0)

  # Load in one file to get the SNPs.
  load(files[1])
  rownames(prsmth) = make.names(rownames(prsmth))
  snps = snps[snps[,1] %in% rownames(prsmth),]

  results = data.frame(SNP = snps[,1], matrix(0, nrow(snps), 8 + 36,
            dimnames = list(snps[,1], c(colnames(prsmth), LETTERS[1:8]))))

  for(f in files) {
    print(f)
    # This loads 'prsmth' that contains genotype probabilities on a log scale.
    load(f)
    results[,2:37] = results[,2:37] + exp(prsmth)
  } # for(f)
  results[,2:37] = results[,2:37] / rowSums(results[,2:37])

  # Summarize founder probabilities.
  mat = matrix(0, ncol(prsmth), 8, dimnames = list(colnames(prsmth),
               LETTERS[1:8]))
  for(i in 1:ncol(mat)) {
    mat[grep(colnames(mat)[i], rownames(mat)),i] = 0.5
  } # for(i)
  mat[cbind(which(rownames(mat) %in% paste(colnames(mat),
      colnames(mat), sep = "")), 1:ncol(mat))] = 1
  results[,38:45] = as.matrix(results[,2:37]) %*% mat

  return(results)
} # summarize.by.snps()


################################################################################
# Go through the genotype probability files and summarize by samples.
# Arguments: path: character, full path to the directory with the
#                  *.genotype.probs.Rdata files.
#            snps: data.frame with 4 columns containing the SNP ID, Chr, Mb &
#                  cM position of each SNP on the genotyping array.
# Returns: data.frame with one row per SNP and one column for each genotype
#          and founder.
summarize.by.samples = function(path = ".", snps) {

  path = add.slash(path)
  files = dir(path = path, pattern = "genotype.probs.Rdata", full.names = T)

  if(length(files) == 0) {
    stop(paste("There are no *.genotype.probs.Rdata files in this directory."))
  } # if(length(files) == 0)

  # Load in one file to get the SNPs.
  load(files[1])
  snps = snps[snps[,1] %in% rownames(prsmth),]

  results = data.frame(sample = gsub(paste(path, "|\\.genotype.probs.Rdata",
            sep = ""), "", files), matrix(0, length(files), 8 + 36,
            dimnames = list(NULL, c(colnames(prsmth), LETTERS[1:8]))))
  index = 1
  for(f in files) {
    print(f)
    # This loads 'prsmth' that contains genotype probabilities on a log scale.
    load(f)
    prsmth = exp(prsmth)
    results[index,2:37] = colMeans(prsmth)
    index = index + 1
  } # for(f)

  # Summarize founder probabilities.
  mat = matrix(0, ncol(prsmth), 8, dimnames = list(colnames(prsmth),
               LETTERS[1:8]))
  for(i in 1:ncol(mat)) {
    mat[grep(colnames(mat)[i], rownames(mat)),i] = 0.5
  } # for(i)
  mat[cbind(which(rownames(mat) %in% paste(colnames(mat),
      colnames(mat), sep = "")), 1:ncol(mat))] = 1
  results[,38:45] = as.matrix(results[,2:37]) %*% mat

  return(results)
} # summarize.by.samples()


################################################################################
# Plot the genotype probabilities across all SNPs.
#
genotype.by.snp.barplot = function(snp.results) {

  data(snps)
  data(founder.cols)

  snps = snps[snps[,1] %in% snp.results$SNP,]
  chr.bnd = table(snps[,2])
  chr.bnd = cumsum(chr.bnd[order(as.numeric(names(chr.bnd)))])
  chr.mid = c(0, chr.bnd[-length(chr.bnd)]) + diff(c(0, chr.bnd)) * 0.5

  layout(matrix(1:8, 8, 1))
  par(font = 2, font.lab = 2, font.axis = 2, las = 1,
      plt = c(0.01, 0.95, 0.05, 0.95))
  ylim = c(0, max(snp.results[,38:45]))
  for(i in 1:8) {
    barplot(snp.results[,37+i], col = founder.cols[i,3], border = 
            founder.cols[i,3], space = 0, ylim = ylim, axes = F)
    axis(2, at = 0:3 * 0.1, pos = 0, cex = 1.25)
    abline(v = chr.bnd, col = "grey80")
    mtext(founder.cols[i,2], side = 4, line = -1.5)
    text(chr.mid, ylim[2] * 0.95, names(chr.bnd), cex = 1.5)
  } # for(i)

} # genotype.by.snp.barplot()


################################################################################
# Plot the genotype probabilities across all samples.
#
genotype.by.sample.barplot = function(results) {

  data(founder.cols)

  # Set up the plot, making breakpoints and converting the founder matrix to 
  # a difference matrix.
  par(font = 2, font.lab = 2, font.axis = 2, las = 1,
      plt = c(0.05, 0.95, 0.05, 0.95))
  mat = t(as.matrix(results[,38:45]))
  colnames(mat) = results$sample
  mat = mat - 0.125
  breaks = quantile(mat, 0:1000 / 1000)
  col = c(colorRampPalette(c(rgb(0,0,1), rgb(0,0,0)))(500),
          colorRampPalette(c(rgb(0,0,0), rgb(1,0,0)))(500))

  image(1:nrow(mat), 1:ncol(mat), mat, ann = F, axes = F, breaks = breaks,
        col = col)
  
  # Write the sample names.
  par(xpd = NA)
  text(0.5, 1:ncol(mat), colnames(mat), adj = 1)
  text(1:nrow(mat), 0, rownames(mat), adj = c(0.5, 1))
} # genotype.by.sample.barplot()


################################################################################
# Make a boxplot of the number of recombinations per mouse by generation.
# You must supply the generation of each sample.
plot.num.recomb = function(results, gen) {

  samples = unique(results$sample)

  if(missing(gen)) {
    gen = rep(1, length(samples))
    names(gen) = samples
  } else {
    if(is.null(names(gen))) {
      stop("gen must have names to line up the sample IDs")
    } # 
  } # else

  tbl = table(results$sample)
  tbl = tbl[names(tbl) %in% names(gen)]
  gen = gen[names(gen) %in% names(tbl)]
  gen = gen[match(names(tbl), names(gen))]

} # plot.num.recomb()



