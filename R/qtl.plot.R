################################################################################
# contains plot.doqtl, effect.plot, pxg.plot
#
# Given a set of mapping statistics and a set of SNPs, make a simple QTL plot.
# Arguments: x: list containing two elements from scanone.
#                   lod: data.frame with chr locations and LOD score.
#                   coef: matrix with model coeffieints.
#            stat.name: character, the name of the mapping statistic to plot.
#                       Must be one of "lod" or "neg.log10.p".
#            sig.thr: matrix, numeric, a set of significance thresholds to plot.
#                     Columns should be names "A" and "X" in that order and
#                     should contain threhsolds for each chromosome obtained
#                     from get.sig.thr().
#            sig.col: vector, color, a set of colors to use for each
#                     significance threshold. Must be same length as sig.thr.
# To plot a subset of chromosomes, feed in a subset of SNPs.
plot.doqtl = function(x, stat.name = c("lod", "neg.log10.p"),  sig.thr = NULL, 
             sig.col = "red", ...) {
  
  doqtl = x
  chrlen = get.chr.lengths()
  stat.name = match.arg(stat.name)
  # Get the call and arguments.
  call = match.call()
  lod = doqtl$lod$A
  if(any(names(doqtl$lod) == "X")) {
    lod = rbind(doqtl$lod$A, doqtl$lod$X)
  } # if(length(lod) > 1)
  # Get chr lengths and locations.  Create an x-axis based on Genome Mb.
  if(max(lod[,3], na.rm = TRUE) > 200) {
    lod[,3] = lod[,3] * 1e-6
  } # if(max(lod[,3]) > 200)
  mb = lod[,3]
  gmb = mb
  unique.chr = as.character(unique(lod[,2]))
  old.warn = options("warn")$warn
  options(warn = -1)
  unique.chr = unique.chr[order(as.numeric(unique.chr))]
  options(warn = old.warn)
  chrlen = chrlen[names(chrlen) %in% unique.chr]
  if(max(chrlen) > 200) {
    chrlen = chrlen * 1e-6
  } # if(max(chrlen) > 200)
  # Get the cumulative sum of the Chr lengths.
  chrsum = cumsum(chrlen)
  chrsum = c(0, chrsum)
  chrmid = chrsum[-length(chrsum)] + (diff(chrsum) * 0.5)
  # Add the preceeding chromosome lengths to each SNP position.
  for(c in 2:length(unique.chr)) {
    rows = which(lod[,2] == unique.chr[c])
    gmb[rows] = gmb[rows] + chrsum[c]
  } # for(c)
  # Make the basic plot.
  plot.column = which(colnames(lod) == stat.name)
  if(length(plot.column) == 0) {
    stop(paste("The stat.name of", stat.name, "was not found in the column names",
         "of the DOQTL object. Please verify that stat.name contains one of the",
         "column names in the DOQTL object."))
  } # if(length(plot.column) == 0)
  par(font = 2, font.lab = 2, font.axis = 2, las = 1, xaxs = "i",
      plt = c(0.12, 0.95, 0.05, 0.88))
  if("ylim" %in% names(call)) {
    plot(gmb, lod[,plot.column], col = 0, xlab = "", xaxt = "n", ylab = stat.name, ...)
  } else {
    if(!missing(sig.thr)) {
      plot(gmb, lod[,plot.column], col = 0, xlab = "", ylab = stat.name,
           xaxt = "n", ylim = c(0, max(lod[,plot.column], sig.thr, na.rm = TRUE) * 1.05), ...)
    } else {
      plot(gmb, lod[,plot.column], col = 0, xlab = "", ylab = stat.name,
           xaxt = "n", ylim = c(0, max(lod[,plot.column], na.rm = TRUE) * 1.05), ...)
    } # else
  } # else
  lod = cbind(lod, gmb)
  lod = split(lod, lod[,2])
  usr = par("usr")
  rect(chrsum[2 * 1:(length(chrsum) * 0.5) - 1], usr[3],
       chrsum[2 * 1:(length(chrsum) * 0.5)], usr[4], col = rgb(0.8,0.8,0.8),
       border = NA)
  rect(usr[1], usr[3], usr[2], usr[4], border = 1)
  text(chrmid, 0.95 * usr[4], names(chrsum)[-1])
  lapply(lod, function(z) { points(z$gmb, z[,plot.column], type = "l", lwd = 2)})
  if(!is.null(sig.thr)) {
    add.sig.thr(sig.thr = sig.thr, sig.col = sig.col, chrsum = chrsum)
  } # if(!is.null(sig.thr))
} # plot.doqtl()
################################################################################
# Using the 8 state model data, make and effect plot that shows the coefficients
# of each founder strain and the LOD.
# Daniel Gatti
# Dan.Gatti@jax.org
# June 1, 2011
# Aug. 4, 2013: Changed to accept doqtl object with separate coefficients on 
#               X chr.
################################################################################
# Arguments: doqtl: a list, as output from scanone, containing two elements. 
#                   lod: data.frame with 7 columns containing the chromosome
#                        locations and lod score.
#                   coef: numeric matrix containing the model coefficients at 
#                         the same loci as the qtl.  Rownames should be the 
#                         same as lod[,1] and column names should be named for 
#                         the founders.
#            chr: character, the chromosome to plot.
#            stat.name: character, name of the mapping statistic.
#            conf.int: boolean that is true if the confidence interval should
#                      be shaded.
#            legend: boolean, true if legend should be drawn.
#            colors: if plotting DO or CC mice, use "DO". This will place the
#                    standard DO colors. If not, a data.frame containing the
#                    founder letter in column 1, the founder name in column 2,
#                    and the colors to use in column 3. See data(founder.colors)
#                    for an example of the DO colors.
#            sex: characters, either "M" or "F". Only used when plotting the X 
#                 chromosome.
coefplot = function(doqtl, chr = 1, stat.name = "LOD", conf.int = TRUE, legend = TRUE,
                colors = "DO", sex, ...) {
  old.par = par(no.readonly = TRUE)
  cross = attr(doqtl, "cross")
  if(is.null(cross)) {
    if(colors[1] == "DO") {    
      colors = do.colors
    } else if(colors[1] == "HS") {
      colors = hs.colors
    } # else if(colors[1] == "HS")
  } else {
    if(cross == "DO") {    
      colors = do.colors
    } else if(cross == "HS") {
      colors = hs.colors
    } # else if(cross == "HS")
  } # else
  num.founders = nrow(colors)
  call = match.call()
  # Keep only the founder coefficients from the coef.matrix.
  lod = NULL
  coef = NULL
  if(chr == "X") {
    if(missing(sex)) {
      stop("Sex (either M or F) must be specified on X chromosome.")
    } # if(missing(sex))
    lod  = doqtl$lod$X
    coef = doqtl$coef$X
    if(sex == "F") {
      columns = match(paste("F", colors[,1], sep = "."), colnames(coef))
    } else {
      columns = match(paste("M", colors[,1], sep = "."), colnames(coef))
    } # else
    columns = columns[!is.na(columns)]
    coef = coef[,c(1,columns)]
    colnames(coef)[1] = "A"
    colnames(coef) = sub("^[MF]\\.", "", colnames(coef))
    coef = coef[rownames(coef) %in% lod[,1],]
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } else {
    lod = doqtl$lod$A
    lod = lod[lod[,2] == chr,]
    intercept = doqtl$coef$A[,1]
    coef = doqtl$coef$A[,(ncol(doqtl$coef$A)-num.founders+1):ncol(doqtl$coef$A)]
    coef[,1] = intercept
    colnames(coef)[1] = "A"
    coef = coef[rownames(coef) %in% lod[,1],]
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } # else 
  # Verify that the SNP IDs in the lod & coef matrices match.
  if(!all(lod[,1] == rownames(coef))) {
    stop(paste("The SNP IDs in column 1 of the qtl data frame must match",
         "the SNP IDs in the rownames of the coef matrix."))
  } # if(!all(lod[,1] == rownames(coef)))
  # Verify that the coefficient column names are in the colors.
  if(!all(colnames(coef) %in% colors[,1])) {
    stop(paste("The founder names in the colnames of the coefficient matrix",
         "must be in column 1 of the colors matrix."))
  } # if(!all(colnames(coef) %in% colors[,1]))
  # Convert the chromosome locations to Mb.
  if(max(lod[,3], na.rm = TRUE) > 200) {
    lod[,3] = lod[,3] * 1e-6
  } # if(max(lod[,3], na.rm = TRUE) > 200)
  # Set the layout to plot the coefficients on top and the p-values on the 
  # bottom.
  layout(mat = matrix(1:2, 2, 1), heights = c(0.66666667, 0.3333333))
  par(font = 2, font.lab = 2, font.axis = 2, las = 1, plt =
      c(0.12, 0.99, 0, 0.85), xaxs = "i", lwd = 2)
  # Plot the coefficients.
  plot(lod[,3], coef[,colors[1,1]], type = "l", col = colors[1,3], lwd = 2,
       ylim = c(min(coef, na.rm = TRUE), max(coef * 2, na.rm = TRUE)), xlab = 
       paste("Chr", chr), ylab = "Founder Effects", axes = FALSE, ...)
  abline(v = 0:20 * 10, col = "grey80")
  for(i in 1:nrow(colors)) {
    points(lod[,3], coef[,colors[i,1]], type = "l", col = colors[i,3],
           lwd = 2)
  } # for(i)
  # Draw a legend for the founder names and colors.
  if(legend) {
    legend.side = "topleft"
    if(which.max(lod[,7]) < nrow(lod) / 2) {
      legend.side = "topright"
    } # if(which.max(apply(coef, 1, max)) < nrow(lod) / 2)
    legend(legend.side, colors[,2], col = colors[,3], lty = 1, lwd = 2,
           x.intersp = 0.75, y.intersp = 0.75, bg = "white", cex = 0.8)
  } # if(legend)
  # Add the axis.
  axis(2)
  # Plot a rectangle around the plot.
  par(xpd = NA)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(xpd = FALSE)
  # Plot the mapping statistic.
  par(plt = c(0.12, 0.99, 0.35, 1))
  # Single chromosome plot.
  plot(lod[,3], lod[,7], type = "l", lwd = 2, xlab = "",
       ylab = stat.name, ...)
  abline(v = 0:20 * 10, col = "grey80")
  points(lod[,3], lod[,7], type = "l", lwd = 2)
  # Shade the confidence interval.
  if(conf.int) {
    interval = bayesint(doqtl, chr = chr)
    usr = par("usr")
    rect(interval[1,3], usr[3], interval[3,3], usr[4], col = rgb(0,0,1,0.1), 
         border = NA)
  } # if(!is.na(conf.int))
  mtext(paste("Chr", chr), 1, 2)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(old.par)
} # coefplot()
# Plot the 36 state effect plot.
# Arguments: pheno: data.frame with with phenotype data.  Samples in rows, 
#                   phenotypes in columns.
#            col: numeric index into phenotype matrix of the phenotype to plot.
#            probs: 3D numeric array, with founder probabilities for
#                  all SNPs and samples. Samples in dim[1], founders in dim[2],
#                  SNPs in dim[3].
#            snp.id: character or numeric, the SNP ID of the SNP to use in the
#                 plot. If numeric, this is the index of the SNP to plot.
#                 If character this is the SNP name to plot.
#            legend: boolean, if TRUE, then add a legend. Default = TRUE.
#            sex: character vector, with the sex of each sample as M or F.
#                 This is used only when the SNP ID is on the X chromosome
#                 and all of the samples are male.  In this case, there are 
#                 only 8 genotype states.
#            covar: a vector that can be converted to a factor with 
#                   a category (i.e. sex or diet, etc.) for each sample. Must
#                   be named with the sample IDs that match rownames(pheno).
#            ...: arguments passed along to plot.
pxg.plot = function(pheno, pheno.col, probs, snp.id, snps, legend = TRUE,
               sex = NA, covar, ...) {
  if(is.null(pheno)) {
    stop(paste("The phenotype matrix cannot be null."))
  } # if(is.null(pheno))
  if(is.null(rownames(pheno))) {
    stop(paste("The phenotypes must have rownames to verify that samples",
         "are lined up."))
  } # if(is.null(rownames(pheno)))
  if(is.null(pheno.col)) {
    stop(paste("The phenotype column cannot be null."))
  } # if(is.null(pheno.col))
  if(is.null(probs)) {
    stop(paste("The founder probabilities array cannot be null."))
  } # if(is.null(probs))
  pheno = pheno[rownames(pheno) %in% dimnames(probs)[[1]],,drop=FALSE]
  probs = probs[dimnames(probs)[[1]] %in% rownames(pheno),,]
  pheno = pheno[match(dimnames(probs)[[1]], rownames(pheno)),,drop=FALSE]
  if(nrow(pheno) == 0 | dim(probs)[1] == 0) {
    stop(paste("The sample names in pheno and probs do not match."))
  } # if(nrow(pheno) == 0 | dim(probs)[1] == 0)
  if(any(rownames(pheno) != dimnames(probs)[[1]])) {
    stop(paste("The phenotypes must have the same sample names as",
         "the founder probabilities."))
  } # if(nrow(pheno) != dim(probs)[[1]])
  if(!missing(covar)) {
    covar = covar[rownames(pheno)]
    if(!all(names(covar) %in% rownames(pheno))) {
      stop(paste("The names of covar must be equal to the rownames of pheno."))
    } # if(all(names(covar), rownames(pheno)))
    covar = factor(covar)
  } # if(!missing(covar))
  if(is.null(snp.id)) {
    stop(paste("The snp.id cannot be null."))
  } # if(is.null(snp.id))
  if(is.null(snps)) {
    stop(paste("The snps cannot be null."))
  } # if(is.null(snps))
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  rownames(snps) = snps[,1]
  # Get the phenotype and genotype probabilities.
  slice = probs[,,snp.id]
  ph    = pheno[,pheno.col]
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
      # If one genotype is > 0.75, then sample is homozygous.
      if(max(slice[i,]) > 0.75) {
        mx = which.max(slice[i,])
        gt[i] = paste(names(mx)[1], names(mx)[1], sep = "")
      } else {
        mx = colnames(slice)[rank(slice[i,]) > 6]
        gt[i] = paste(sort(mx), collapse = "")
      } # else
    } # for(i)
  } else {
    for(i in 1:nrow(slice)) {
      # If one genotype is > 0.75, then sample is homozygous.
      if(max(slice[i,]) > 0.75) {
        mx = which.max(slice[i,])
        gt[i] = paste(names(mx)[1], names(mx)[1], sep = "")
      } else {
        mx = colnames(slice)[rank(slice[i,]) > 6]
        gt[i] = paste(sort(mx), collapse = "")
      } # else
    } # for(i)
  } # else
  # Order the levels and plot.
  founders = sort(dimnames(probs)[[2]])
  states = outer(founders, founders, paste, sep = "")
  states = sort(states[upper.tri(states, diag = TRUE)])
  # Get means and standard errors.  
  spl = split(ph, gt)
  geno.means = sapply(spl, mean, na.rm = TRUE)
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
  gt = factor(gt, levels = names(geno.means), ordered = TRUE)
  # If we have a covariate, then set the pch to plot the different
  # categories.
  pch = 16
  if(!missing(covar)) {
    pch = as.numeric(covar)
  } # if(!missing(covar))
  # Make the plot.
  par(font = 2, font.lab = 2, font.axis = 2, las = 2)
  # If we're plotting the legend, add a little space at the top of the plot.
  ylim = range(ph, na.rm = TRUE)
  if(legend) {
    ylim[2] = ylim[2] + (ylim[2] - ylim[1]) * 0.2
  } # if(legend)
  plot(1, 1, col = 0, xlim = c(0.5, length(states) + 0.5),
       ylim = ylim, xaxt = "n", xlab = "Genotype",
       ylab = colnames(pheno)[pheno.col], ...)
  abline(v = 1:36, col = "grey80")
  points(as.numeric(gt), ph, pch = pch, col = rgb(0,0,0,0.4))
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
    legend("topleft", pch = do.colors[,1], legend = do.colors[,2],
           bg = "white", cex = 0.9)
  } # if(legend)
} # pxg.plot()
