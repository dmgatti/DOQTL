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
coef.plot = function(doqtl, chr = 1, stat.name = "LOD", conf.int = T, legend = T,
            colors = "DO", sex, ...) {

  old.par = par(no.readonly = TRUE)

  if(colors[1] == "DO") {
    data(founder.cols)
    colors = founder.cols
  } # if(colors == "DO")

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
  if(max(lod[,3], na.rm = T) > 200) {
    lod[,3] = lod[,3] * 1e-6
  } # if(max(lod[,3], na.rm = T) > 200)

  # Set the layout to plot the coefficients on top and the p-values on the 
  # bottom.
  layout(mat = matrix(1:2, 2, 1), heights = c(0.66666667, 0.3333333))
  par(font = 2, font.lab = 2, font.axis = 2, las = 1, plt =
      c(0.12, 0.99, 0, 0.85), xaxs = "i", lwd = 2)

  # Plot the coefficients.
  plot(lod[,3], coef[,colors[1,1]], type = "l", col = colors[1,3], lwd = 2,
       ylim = c(min(coef, na.rm = T), max(coef * 2, na.rm = T)), xlab = 
       paste("Chr", chr), ylab = "Coefficient", axes = F, ...)
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
  par(xpd = F)

  # Plot the mapping statistic.
  par(plt = c(0.12, 0.99, 0.3, 1))
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
} # coef.plot()

