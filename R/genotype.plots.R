################################################################################
# Make a plot of DO genotypes on each chromosome.
# Feb. 14, 2011
# Daniel Gatti
# Dan.Gatti@jax.org
# Contains: chr.skeletons, genomic.points, plot.genoprob,
#           write.genoprob.plots.
################################################################################

#############################################################################
# Draw the skeleton of the chromosomes.
# Arguments: chr: vector with chr names to plot.
chr.skeletons = function(chr, chrlen, ...) {

  # Plot a coordinate system to draw on.
  par(font = 2, font.axis = 2, font.lab = 2, las = 1, plt = c(0.14, 0.95, 
      0.1, 0.9))
  plot(1, 1, col = 0, xlim = c(1, (2 * length(chrlen) + 1)),
       ylim = c(0, 200), xaxt = "n", yaxs = "i", xlab = "",
       ylab = "Mb")

  # Draw chr skeletons
  for(i in 1:length(chrlen)) {
    lines(2 * c(i, i), c(0, chrlen[i]))
  } # for(i)

  # Draw light lines at 10 Mb and heavier lines at 50 Mb.
  abline(h = 0:20 * 10, col = "grey80")
  abline(h = 0:4  * 50, col = "grey50", lwd = 1.2)
  mtext(side = 1, at = 2 * 1:length(chrlen), text = names(chrlen))

} # chr.skeletons()


################################################################################
# Once a Chr skeleton is up, add points to it.
# Arguments: chr: numeric vector, with chr for each point to plot.
#            loc: numeric vector, with Mb locations for each point to plot.
#                 Must be the same length as chr.
genomic.points = function(chr = NULL, loc = NULL) {
  if(!is.null(chr) & !is.null(loc)) {
    if(length(chr) != length(loc)) {
      stop(paste("chr and loc are not the same length. Please use chr and",
                 "loc vectors of the same length."))
    } # if(length(chr) != length(loc))
    points(2 * chr, loc, pch = 3, cex = 2, lwd = 2, col = 2)
  } # if(!is.null(chr) & !is.null(loc))
} # genomic.points()


################################################################################
# Plot the genotypes using only the maximum probability at each SNP and arrange
# the colors to minimize jumping from one strand to the other.
# Arguments: x: numeric matrix, with smoothed probabilities. SNPs in rows, 
#               genotypes in columns. May also be a genoprobs object with
#               a cross attribute specifying DO or HS mice.
#            snps: matrix with SNP IDs in column 1, chromosomes in column 2 and
#                  Mb locations in column 3.
#            colors: data.frame with three columns containing the founder 
#                    letters, names and colors in columns 1, 2 & 3, respectively.
#                    Default = "DO", which automatically uses the CC founder 
#                    colors.
#            chrlen: named numeric vector containing chromosome lengths or "mm10" for the mouse chromosomes.
#            ...: other arguments to be passed to plot.
plot.genoprobs = function(x, snps, colors = "DO", chrlen = "mm10", ...) {

  old.par = par(no.readonly = TRUE)

  states = colnames(x)
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3], na.rm = TRUE) > 200)
  
  cross = attr(x, "cross")
  if(is.null(cross)) {
    if(colors == "DO") {
      colors = do.colors
    } else {
      if(!is.data.frame(colors)) {
        stop(paste("If not using the DO colors, colors must be a data frame",
  	         "with three columns containing the founder letters, names and",
  		       "colors in columns 1, 2 & 3, respectively."))
      } # if(!is.data.frame(colors))
	
      if(ncol(colors) != 3) {
        stop(paste("Colors does not have three columns. Colors must have three",
	         "columns containing the founder letters, names and colors in",
          "columns 1, 2 & 3, respectively."))
	} # if(ncol(colors) != 3)
    } # else
  } else {
    if(cross == "DO" | cross == "CC" | cross == "DOF1") {
      colors = do.colors
    } else if(cross == "HS") {
      colors = hs.colors 
    } # else
  } # else

  if(chrlen == "mm10") {
    chrlen = get.chr.lengths(chrlen)
  } # if(chrlen == "DO")

  # Subset the data to only include the SNPs in the snpfile.
  x = x[rownames(x) %in% snps[,1],]
  snps = snps[snps[,1] %in% rownames(x),]
  prmsth = x[match(snps[,1], rownames(x)),]

  # Create state colors.
  state.cols = cbind(states, matrix(unlist(strsplit(states, split = "")),
               length(states), 2, byrow = TRUE, dimnames = list(states, NULL)))
  state.cols[,2] = colors[match(state.cols[,2], colors[,1]),3]
  state.cols[,3] = colors[match(state.cols[,3], colors[,1]),3]
  # If we see an "I" in the states, then this is DO/F1 data and we color
  # the other founder in dark grey.
  if(length(grep("I", state.cols[,1])) > 0) {
    state.cols[,3] = rgb(0.3, 0.3, 0.3)
  } # if(length(grep("I", state.cols[,1])) > 0)
  # Plot the chromosome skeletons.
  old.warn = options("warn")$warn
  options(warn = -1)
  chr = sort(as.numeric(unique(snps[,2])), na.last = TRUE)
  options(warn = old.warn)
  chr[20] = "X"
  chrlen = get.chr.lengths()
  par(plt = c(0.08, 0.95, 0.05, 0.88))
  chr.skeletons(chr, chrlen, ...)
  chrlen = chrlen[names(chrlen) %in% snps[,2]]
  
  # Get the maximum probabilities.
  maxstate = state.cols[apply(x, 1, which.max),1]
  names(maxstate) = rownames(x)
  offset = 0.6
  for(c in 1:length(chrlen)) {
    ss = which(snps[,2] == names(chrlen)[c])
    if(length(ss) > 0) {
      ms = maxstate[ss]
      # xm is the middle of the chromosome.
      xm = 2 * c
      # xl is the left side.
      xl = xm - offset
      # xr is the right side.
      xr = xm + offset
      y = snps[ss,3]
      # Get the left and right colors.
      lcol = state.cols[ms,2]
      rcol = state.cols[ms,3]
      # Get the left and right side breaks points.
      lbreaks = c(1, which(lcol[1:(length(lcol)-1)] != lcol[2:length(lcol)]) +
                  1, length(lcol))
      rbreaks = c(1, which(rcol[1:(length(rcol)-1)] != rcol[2:length(rcol)]) +
                  1, length(rcol))
      lcol = lcol[lbreaks]
      rcol = rcol[rbreaks]
      # Make a data.frame of all breaks and try to minimize swapping colors
      # from one chr to the other.
      all.breaks = sort(union(lbreaks, rbreaks))
      all.breaks = data.frame(breaks = all.breaks, left = rep(NA,
                   length(all.breaks)), right = rep(NA, length(all.breaks)))
      all.breaks$left[match(lbreaks, all.breaks$breaks)]  = lcol
      all.breaks$right[match(rbreaks, all.breaks$breaks)] = rcol
      # Fill in any NA values.
      for(i in 1:nrow(all.breaks)) {
        if(is.na(all.breaks$left[i])) {
          all.breaks$left[i] = all.breaks$left[i-1]
        } # if(is.na(all.breaks$left[i]))
        if(is.na(all.breaks$right[i])) {
          all.breaks$right[i] = all.breaks$right[i-1]
        } # if(is.na(all.breaks$right[i])) 
      } #for(i)
      rect(xl, y[all.breaks$breaks[1]], xm, y[all.breaks$breaks[2]],
          col = all.breaks$left[1], density = NA, border =
          all.breaks$left[1])
      # Plot the right side of the chromosome.
      rect(xm, y[all.breaks$breaks[1]], xr, y[all.breaks$breaks[2]],
           col = all.breaks$right[1], density = NA, border =
           all.breaks$right[1])
      # Go through each break and see if the next color should be swapped.
      for(i in 2:nrow(all.breaks)) {
        if(all.breaks$left[i] != all.breaks$right[i]) {
          if(all.breaks$left[i] == all.breaks$right[i-1] |
             all.breaks$left[i-1] == all.breaks$right[i]) {
               tmp = all.breaks$left[i]
               all.breaks$left[i] = all.breaks$right[i]
               all.breaks$right[i] = tmp
          } # if(all.breaks$left[i] == all.breaks$right[i-1] | ...
        } # if(all.breaks$left[i] != all.breaks$right[i])
      	
   	  rect(xl, y[all.breaks$breaks[i]], xm, y[all.breaks$breaks[i+1]],
   	     col = all.breaks$left[i], density = NA, border =
   	     all.breaks$left[i])
        # Plot the right side of the chromosome.
        rect(xm, y[all.breaks$breaks[i]], xr, y[all.breaks$breaks[i+1]],
        col = all.breaks$right[i], density = NA, border =
    	  all.breaks$right[i])
      } # for(i)
    } # if(length(ss) > 0)
  } # for(c)
  legend(x = 28, y = 200, legend = colors[,2], col = colors[,3], pch = 15,
         bg = "white")
  
} # plot.genoprobs()
################################################################################
# Make genotype plots for all of the files in a directory.
# This requires that the *.Rdata files have been made.
write.genoprob.plots = function(path = ".", snps) {
  files = dir(path, pattern = "genotype.probs.Rdata", full.names = TRUE)
  if(!is.null(files)) {
    prsmth = NULL
    for(f in files) {
      print(f)
      # This load in x.
      load(f)
      sample = gsub(paste("^", path, "/|\\.genotype\\.probs\\.Rdata$", sep = ""),
                    "", f)
      png(paste(sample, "genotype.probs.png", sep = "."), width = 1000,
          height = 1000, res = 144)
      plot.genoprobs(x = prsmth, snps, main = sample)
      dev.off()
    } # for(f)
  } # if(!is.null(files))
} # write.genoprob.plots()
