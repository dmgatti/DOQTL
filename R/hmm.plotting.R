################################################################################
# Functions for plotting of DO genotypes.
# Jan. 11, 2012
# Daniel Gatti
# Dan.Gatti@jax.org
################################################################################
options(stringsAsFactors = FALSE)
################################################################################
# Plot the smoothed probability matrix for the requested sample.
# Arguments: sample.index: integer, index into prsmth array of sample to plot.
#            states: character vector, DO genotype state names.
#            prsmth: numeric 3D array, smoothed probabilities for all samples.
prsmth.plot = function(sample.index, states, prsmth) {
  par(yaxt = "n", font = 2, font.axis = 2, font.lab = 2, las = 1)
  breaks = 0:101/100
  col = rev(grey(breaks[-length(breaks)]))
  image(1:dim(prsmth)[3], 1:dim(prsmth)[1],
        exp(t(prsmth[dim(prsmth)[1]:1,sample.index,])),
        xlab = "SNPs", ylab = "Genotypes", main =
        dimnames(prsmth)[[2]][sample.index], breaks = breaks, col = col)
  abline(h = 0:length(states) + 0.5, col = "grey70")
  abline(h = c(0,length(states)) + 0.5, col = 1)
  abline(v = c(0, dim(prsmth)[[3]]) + 0.5, col = 1)
  par(las = 2)
  mtext(text = rev(states), side = 2, line = 0.2, at = 1:length(states))
  mtext(text = rev(states), side = 4, line = 0.2, at = 1:length(states))
} # prsmth.plot()
################################################################################
# Plot the means and variances for the requested SNP.
# Arguments: s: integer, the index into the SNPs..
#            states: vector, chracter, the HMM states.
#            theta: the thete intensity matrix.
#            rho: the rho intensity matrix.
#            r.t.means: the state mean matrix for the current chromosome.
#            r.t.covars: the state covariance matrices for the current
#                        chromosome.
#            sample: (optional) integer that is the index of a sample to plot.
intensity.mean.covar.plot = function(s, states, theta, rho, r.t.means, r.t.covars,
                            sample, ...) {
  
  par(font = 2, font.axis = 2, font.lab = 2, yaxt = "s", las = 1)
  col = do.colors
  xlim = range(theta[,s], r.t.means[,1,s], na.rm = TRUE)
  ylim = range(rho[,s],   r.t.means[,2,s], na.rm = TRUE)
  plot(theta[,s], rho[,s], pch = 16, main = paste(s, colnames(theta)[s]),
       xlab = expression(theta), ylab = expression(rho), xlim = xlim, 
       ylim = ylim, ...)
  
  # If the states passed in are hemizygous, change them to 2 letters.
  if(nchar(states[1]) == 1) {
    states = paste(states, states, sep = "")
  } # if(nchar(states[1]) = 1)
  state.cols = matrix(unlist(strsplit(states, split = "")), nrow = 2)
  for(i in 1:nrow(col)) {
    state.cols[state.cols == col[i,1]] = col[i,3]
  } # for(i)
  # Plot the covariances.
  for(i in 1:dim(r.t.covars)[1]) {
    lines(c(r.t.means[i,1,s] - sqrt(r.t.covars[i,1,s]),
            r.t.means[i,1,s] + sqrt(r.t.covars[i,1,s])),
            c(r.t.means[i,2,s], r.t.means[i,2,s]), col = 
            "light grey")
    lines(c(r.t.means[i,1,s], r.t.means[i,1,s]),
          c(r.t.means[i,2,s] - sqrt(r.t.covars[i,2,s]), 
            r.t.means[i,2,s] + sqrt(r.t.covars[i,2,s])),
            col = "light grey")
  } # for(i)
  
  # Plot the means.
  points(r.t.means[,1,s], r.t.means[,2,s], pch = 16, col =
         state.cols[1,], cex = 2.2)
  points(r.t.means[,1,s], r.t.means[,2,s], pch = 17, col =
         state.cols[2,], cex = 1.5)
  # If provided, plot the sample.
  if(!missing(sample)) {
    points(theta[sample,s], rho[sample,s], pch = 16, cex = 2, col = "orange")
    points(theta[sample,s], rho[sample,s], pch = 4,  cex = 3, lwd = 3,
           col = "orange")
  } # if(!missing(sample))
} # intensity.mean.covar.plot()
##############################################################################
# Plot the founders and F1s at a given SNP.
# This assumes that the colnames on theta contain some two letter founder IDs.
# Arguments: theta: the matrix of theta intensites with SNPs in rows and samples
#               in columns.  Founder/F1s must be named with 2 letter codes.
#            rho: the matrix of rho intensites with SNPs in rows and samples
#               in columns.  Founder/F1s must be named with 2 letter codes.
#            s: the index of the SNP to plot.
#            is.founder.F1: vector of T/F values that is T if a sample is a
#                           founder or F1.
#            plotNew: T/F that is T if a new plot should be drawn or F if the
#                     founders should be added to an existing plot.
#            ...: Additional graphical parameters passed to plot.
founder.F1.intensity.plot = function(theta, rho, s = 1, is.founder.F1, plotNew = TRUE, ...) {
  if(all(is.founder.F1 == FALSE)) {
    stop("!!! No founders or F1s in plot.founders() !!!")
  } # if(all(is.founder.F1 == FALSE))
  if(plotNew) {
    plot(theta[,s], rho[,s], main = colnames(theta)[s], pch = 16, ...)
  } # if(plotNew)
  do.colors[,1] = paste(do.colors[,1], do.colors[,1], sep = "")
  spl = strsplit(rownames(theta)[is.founder.F1], split = "")
  spl = sapply(spl, function(a) { paste(a, a, sep = "") })
  spl = matrix(do.colors[match(unlist(spl), do.colors[,1]),3], nrow = 2,
        dimnames = list(NULL, rownames(theta)[is.founder.F1]))
  points(theta[is.founder.F1, s], rho[is.founder.F1, s], pch = 16, cex = 2.2,
         col = spl[1,])
  points(theta[is.founder.F1, s], rho[is.founder.F1, s], pch = 17, cex = 1.5,
         col = spl[2,])
} # founder.F1.intensity.plot()
