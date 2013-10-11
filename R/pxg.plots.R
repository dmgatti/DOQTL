################################################################################
# Phenotype by genotype plots, with {ACGT} SNPs, founder probabilities (8 state)
# or genotype probabilities (36 state).
# Daniel Gatti
# Dan.Gatti@jax.org
# Feb. 4, 2012
################################################################################
# Plot the phenotype vs. the {ACGT} genotype at one SNP.
# Arguments: pheno: numeric vector, phenotype values for each sample. Sample
#                   names must be in names and match sample names in geno.
#            geno: character matrix, allele calls in {ACGT} format for all SNPs
#                  and samples.  SNPs in rows and samples in columns.
#            snp: character, the SNP ID of the SNP to use in the plot.
pxg.plot = function(pheno, geno, snp) {

  # Make a factor out of the genotypes.
  g = factor(geno[snp,])
  pxg.internal(x = g, y = pheno)

} # pxg.plot()


# Plot the 36 state effect plot.
# Arguments: pheno: data.frame with with phenotype data.  Samples in rows, 
#                   phenotypes in columns.
#            col: numeric index into phenotype matrix of the phenotype to plot.
#            founder.probs: 3D numeric array, with founder probabilities for
#                  all SNPs and samples. Samples in dim[1], founders in dim[2],
#                  SNPs in dim[3].
#            snp.id: character or numeric, the SNP ID of the SNP to use in the
#                 plot. If numeric, this is the index of the SNP to plot.
#                 If character this is the SNP name to plot.
#            legend: boolean, if TRUE, then add a legend. Default = T.
#            sex: character vector, with the sex of each sample as M or F.
#                 This is used only when the SNP ID is on the X chromosome
#                 and all of the samples are male.  In this case, there are 
#                 only 8 genotype states.
effect.plot = function(pheno, col, founder.probs, snp.id, snps, legend = T,
               sex = NA) {

  if(is.null(pheno)) {
    stop(paste("The phenotype matrix cannot be null."))
  } # if(is.null(pheno))

  if(is.null(rownames(pheno))) {
    stop(paste("The phenotypes must have rownames to verify that samples",
         "are lined up."))
  } # if(is.null(rownames(pheno)))

  if(is.null(col)) {
    stop(paste("The phenotype column cannot be null."))
  } # if(is.null(col))

  if(is.null(founder.probs)) {
    stop(paste("The founder probabilities array cannot be null."))
  } # if(is.null(founder.probs))

  pheno = pheno[pheno[,1] %in% dimnames(founder.probs)[[1]],]
  founder.probs = founder.probs[dimnames(founder.probs)[[1]] %in% pheno[,1],,]
  pheno = pheno[match(dimnames(founder.probs)[[1]], pheno[,1]),]

  if(any(pheno[,1] != dimnames(founder.probs)[[1]])) {
    stop(paste("The phenotype must have the same sample names as",
         "the founder probabilities."))
  } # if(nrow(pheno) != dim(founder.probs)[[1]])

  if(is.null(snp.id)) {
    stop(paste("The snp.id cannot be null."))
  } # if(is.null(snp.id))

  if(is.null(snps)) {
    stop(paste("The snps cannot be null."))
  } # if(is.null(snps))

  snps = snps[snps[,1] %in% dimnames(founder.probs)[[3]],]
  rownames(snps) = snps[,1]

  # Get the phenotype and genotype probabilities.
  slice = founder.probs[,,snp.id]
  ph    = pheno[,col]

  # Get the genotype with the maximum probability for each sample.
  gt = rep(NA, nrow(slice))
  # If we are on the X chromosome and we have sex, then process the males
  # and females separately.
  if(snps[snp.id,2] == "X" & !is.na(sex)) {
    sex = toupper(sex)
    males = which(sex == "M")

    # Males
    for(i in males) {
      # Get the maximum genotype.
      mx = which.max(slice[i,])
      gt[i] = paste(names(mx)[1], names(mx)[1], sep = "")
    } # for(i)

    # Females
    females = which(sex == "F")
    for(i in females) {
      # If one genotype is > 0.6, then sample is homozygous.
      if(max(slice[i,]) > 0.6) {
        mx = which.max(slice[i,])
        gt[i] = paste(names(mx)[1], names(mx)[1], sep = "")
      } else {
        mx = colnames(slice)[rank(slice[i,]) > 6]
        gt[i] = paste(sort(mx), collapse = "")
      } # else
    } # for(i)

  } else {
    for(i in 1:nrow(slice)) {
      # If one genotype is > 0.6, then sample is homozygous.
      if(max(slice[i,]) > 0.6) {
        mx = which.max(slice[i,])
        gt[i] = paste(names(mx)[1], names(mx)[1], sep = "")
      } else {
        mx = colnames(slice)[rank(slice[i,]) > 6]
        gt[i] = paste(sort(mx), collapse = "")
      } # else
    } # for(i)
  } # else

  # Order the levels and plot.
  founders = sort(dimnames(founder.probs)[[2]])
  states = outer(founders, founders, paste, sep = "")
  states = sort(states[upper.tri(states, diag = T)])

  # Get means and standard errors.  
  spl = split(ph, gt)
  geno.means = sapply(spl, mean, na.rm = T)
  geno.n = sapply(spl, length)
  geno.se = sapply(spl, sd) / geno.n
  ord = order(geno.means)
  geno.means = geno.means[ord]
  geno.se    = geno.se[ord]

  # Fill in missing genotypes.
  missing.states = states[!states %in% names(geno.means)]
  if(length(missing.states) > 0) {
    geno.means = c(rep(0, length(missing.states)), geno.means)
    geno.se    = c(rep(0, length(missing.states)), geno.se)
    names(geno.means)[1:length(missing.states)] = missing.states
    names(geno.se)[1:length(missing.states)]    = missing.states
  } # if(length(missing.states) > 0)

  # Factor the genotypes.
  gt = factor(gt, levels = names(geno.means), ordered = T)

  # Make the plot.
  par(font = 2, font.lab = 2, font.axis = 2, las = 2)
  # If we're plotting the legend, add a little space at the top of the plot.
  ylim = range(ph, na.rm = T)
  if(legend) {
    ylim[2] = ylim[2] + (ylim[2] - ylim[1]) * 0.2
  } # if(legend)
  plot(1, 1, col = 0, xlim = c(0.5, length(states) + 0.5),
       ylim = ylim, xaxt = "n", xlab = "Genotype",
       ylab = colnames(pheno)[col])
  abline(v = 1:36, col = "grey80")
  points(as.numeric(gt), ph, pch = 20)
  title(main = paste(snps[snp.id,1], "\nChr", snps[snp.id,2], ":",
        snps[snp.id,3], "Mb"))
  axis(side = 1, at = 1:36, labels = names(geno.means))

  # Plot means and standard errors.
  offset = 0.3
  for(i in 1:length(geno.means)) {
    lines(c(i - offset, i + offset), rep(geno.means[i], 2), lwd = 2, col = 2)
    lines(rep(i, 2), c(geno.means[i] - geno.se[i], geno.means[i] + geno.se[i]),
          lwd = 2, col = 2)
  } # for(i)

  # Legend
  if(legend) {
    data(founder.cols)
    legend("topleft", pch = founder.cols[,1], legend = founder.cols[,2],
           bg = "white", cex = 0.9)
  } # if(legend)
} # effect.plot()


# Plot the phenotype vs. the 8 founder state genotypes.
# Arguments: pheno: numeric vector, phenotype values for each sample. Sample
#                   names must be in names and match sample names in geno.
#            founder.probs: 3D numeric array, with founder probabilities for
#                  all SNPs and samples. Samples in dim[1], founders in dim[2],
#                  SNPs in dim[3].
#            snp: character, the SNP ID of the SNP to use in the plot.
pxg.8state.plot = function(pheno, founder.probs, snp) {

  par(mfrow = c(4, 7), plt = c(0.1, 0.9, 0.1, 0.9), font = 2, font.lab = 2,
      font.axis = 2, las = 1)
  g = founder.probs[,,snp]
  for(i in 1:7) {
    for(j in (i+1):8) {
      names = c(paste(colnames(g)[i], colnames(g)[i], sep = ""),
                paste(colnames(g)[i], colnames(g)[j], sep = ""),
                paste(colnames(g)[j], colnames(g)[j], sep = ""))
      AA = which(g[,i] > 0.9)
      AB = which(g[,i] > 0.4 & g[,j] > 0.4)
      BB = which(g[,j] > 0.9)
      x = factor(c(rep(names[1], length(AA)), rep(names[2], length(AB)),
                 rep(names[3], length(BB))))
      y = pheno[c(AA, AB, BB)]
      pxg.internal(x, y)
    } # for(j)
  } # for(i)

} # pxg.8state.plot()


# Internal function for making a simple allele plot.
pxg.internal = function(x, y) {

  # Plot the phenotype values by genotype.
  par(font = 2, font.lab = 2, font.axis = 2,las = 1)
  plot(jitter(as.numeric(x)), y, xaxt = "n", pch = 16, col = "grey50",
       xlab = "Genotype", ylab = "Phenotype", main = snp)

  # Add genotypes to the X-axis.
  mtext(levels(x), side = 1, line = 1, at = 1:length(levels(x)))

  # Add means and Std. Errors.
  p = split(y, x)
  means = sapply(p, mean, na.rm = T)
  se = sapply(p, sd, na.rm = T) / sqrt(sapply(p, length))
  x = 1:length(levels(x))
  for(i in 1:length(x)) {
    lines(c(x[i] - 0.1, x[i] + 0.1), c(means[i], means[i]), col = 2, lwd = 2)
    lines(c(x[i], x[i]), c(means[i] - se[i], means[i] + se[i]), col = 2, 
          lwd = 2)
    lines(c(x[i] - 0.05, x[i] + 0.05), c(means[i] - se[i], means[i] - se[i]),
          col = 2, lwd = 2)
    lines(c(x[i] - 0.05, x[i] + 0.05), c(means[i] + se[i], means[i] + se[i]),
          col = 2, lwd = 2)
  } # for(i)

  abline(v = 0:5 + 0.5, col = "grey80")
} # pxg.internal()
