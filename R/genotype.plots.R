################################################################################
# Make a plot of DO genotypes on each chromosome.
# Feb. 14, 2011
# Daniel Gatti
# Dan.Gatti@jax.org
################################################################################

#############################################################################
# Get the chromosome lengths.
get.chr.lengths = function() {
  # Get the chromosome lengths.
  chrlen = org.Mm.egCHRLENGTHS
  remove = grep("Un|random", names(chrlen))
  if(length(remove) > 0) {
    chrlen = chrlen[-remove]
  } # if(length(remove) > 0)
  old.warn = options("warn")$warn
  options(warn = -1)
  chrlen = chrlen[order(as.numeric(names(chrlen)))]
  options(warn = old.warn)
  chrlen = chrlen * 1e-6
  return(chrlen)
} # get.chr.lengths()


#############################################################################
# Draw the skeleton of the chromosomes.
# Arguments: chr: vector with chr names to plot.
plot.chr.skeletons = function(chr, ...) {
  
  chrlen = get.chr.lengths()
  chrlen = chrlen[names(chrlen) %in% chr]

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

} # plot.chr.skeletons()


################################################################################
# Once a Chr skeleton is up, add points to it.
# Arguments: chr: numeric vector, with chr for each point to plot.
#            loc: numeric vector, with Mb locations for each point to plot.
#                 Must be the same length as chr.
plot.genomic.points = function(chr = NULL, loc = NULL) {

  if(!is.null(chr) & !is.null(loc)) {
    if(length(chr) != length(loc)) {
      stop(paste("chr and loc are not the same length. Please use chr and",
                 "loc vectors of the same length."))
    } # if(length(chr) != length(loc))

    points(2 * chr, loc, pch = 3, cex = 2, lwd = 2, col = 2)
  } # if(!is.null(chr) & !is.null(loc))

} # plot.genomic.points()


################################################################################
# Plot the genotypes using only the maximum probability at each SNP and arrange
# the colors to minimize jumping from one strand to the other.
# Arguments: prsmth: numeric matrix, with smoothed probabilities. SNPs in rows, 
#                    genotypes in columns.
#            snps: matrix with SNP IDs in column 1, chromosomes in column 2 and
#                  Mb locations in column 3.
#            colors: data.frame with three columns containing the founder 
#                    letters, names and colors in columns 1, 2 & 3, respectively.
#                    Default = "DO", which automatically uses the CC founder 
#                    colors.
#            ...: other arguments to be passed to plot.
plot.genotype.max = function(prsmth, snps, colors = "DO", ...) {

  old.par = par(no.readonly = T)
  states = colnames(prsmth)

  if(max(snps[,3], na.rm = T) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3], na.rm = T) > 200)
  
  if(colors == "DO") {
    data(founder.cols)
    colors = founder.cols
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

  # Subset the data to only include the SNPs in the snpfile.
  prsmth = prsmth[rownames(prsmth) %in% snps[,1],]
  snps = snps[snps[,1] %in% rownames(prsmth),]
  prmsth = prsmth[match(snps[,1], rownames(prsmth)),]

  # Create state colors.
  state.cols = cbind(states, matrix(unlist(strsplit(states, split = "")),
               length(states), 2, byrow = T, dimnames = list(states, NULL)))
  state.cols[,2] = founder.cols[match(state.cols[,2], founder.cols[,1]),3]
  state.cols[,3] = founder.cols[match(state.cols[,3], founder.cols[,1]),3]

  # Plot the chromosome skeletons.
  old.warn = options("warn")$warn
  options(warn = -1)
  chr = sort(as.numeric(unique(snps[,2])), na.last = T)
  options(warn = old.warn)
  chr[20] = "X"
  chrlen = get.chr.lengths()
  par(plt = c(0.08, 0.95, 0.05, 0.88))
  plot.chr.skeletons(chr, ...)
  chrlen = chrlen[names(chrlen) %in% snps[,2]]
  
  # Get the maximum probabilities.
  maxstate = state.cols[apply(prsmth, 1, which.max),1]
  names(maxstate) = rownames(prsmth)

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

  legend(x = 30, y = 200, legend = founder.cols[,2], col = founder.cols[,3], pch = 15)
} # plot.genotype.max()


################################################################################
# Make genotype plots for all of the files in a directory.
# This requires that the *.Rdata files have been made.
write.genotype.plots = function(path = ".", snps) {

  files = dir(path, pattern = "genotype.probs.Rdata", full.names = T)

  if(!is.null(files)) {
    data(founder.cols)
    for(f in files) {
      print(f)
      # This load in prsmth.
      load(f)
      sample = gsub(paste("^", path, "/|\\.genotype\\.probs\\.Rdata$", sep = ""),
                    "", f)
      png(paste(sample, "genotype.probs.png", sep = "."), width = 1000,
          height = 1000, res = 144)
      plot.genotype.max(prsmth, snps, title = sample)
      dev.off()
    } # for(f)
  } # if(!is.null(files))
} # make.all.genotype.plots()


# Create genotype plots using the new Viterbi algorithm genotypes.
# Arguments: viterbi: list, with one chr per list element. Each list
#            elelment is a matrix of two letter genotype IDs from the 36
#            genotype states.
#            sample: numeric or character vector. If numeric, the values are
#                    sample indices to plot. If character, the value are
#                    sample names to plot.
plot.genotypes = function(viterbi, sample = 1) {
  data(states)
  data(founder.cols)

  # If we have character samples, then find the corresponding sample indices.
  if(is.character(sample)) {
    sample.char = sample
    sample = match(sample, rownames(viterbi[[1]]))
    if(any(is.na(sample))) {
      print(paste("All of the samples names are not found in the viterbi paths.",
            sample.char[is.na(sample)]))
    } # if(any(is.na(sample)))
  } # if(is.character(sample))
  num.samples = length(sample)

  # Subset the data to keep only the samples that we're plotting.
  for(i in 1:length(viterbi)) {
    viterbi[[i]] = viterbi[[i]][sample,,drop = F]
  } # for(i)

  # Switch the order of the genotypes to somewhat minimize recombinations.
  # This in not intended to be a phasing algorithm.
  for(i in 1:length(viterbi)) { # 1:num.chr
    for(j in 2:ncol(viterbi[[i]])) { # 2:num.snps.on.chr
      if(viterbi[[i]][,j-1] != viterbi[[i]][,j]) {
        let = strsplit(viterbi[[i]][,c(j-1,j)], split = "")
        if(let[[1]][1] != let[[1]][2] & let[[2]][1] != let[[2]][2]) {
          if(let[[1]][1] == let[[2]][2] || let[[1]][2] == let[[2]][1]) {
            viterbi[[i]][,j] = paste(let[[2]][2], let[[2]][1], sep = "")
          }
        }
      } # if(viterbi[[i]][,j-1] != viterbi[[i]][,j])
    } # for(j)
  } # for(i)

  # Get the chr boundaries
  chr = sort(as.numeric(unique(snps[,2])), na.last = T)
  chr[20] = "X"
  chrlen = get.chr.lengths()
  snp.ids = unlist(lapply(viterbi, function(a) { colnames(a) }))
  snps = snps[snps[,1] %in% snp.ids,]
  chrlen = chrlen[names(chrlen) %in% snps[,2]]

  offset = 0.6

  for(i in 1:num.samples) {
    print(paste("Processing Sample", rownames(viterbi[[1]])[i], "..."))
    plot.chr.skeletons(chr)
    title(rownames(viterbi[[1]])[i])
    for(c in 1:length(chrlen)) {
      let = matrix(unlist(strsplit(viterbi[[c]][i,], split = "")), nrow = 2)
      lbrk = which(let[1,-1] != let[1,-ncol(let)])
      rbrk = which(let[2,-1] != let[2,-ncol(let)])
      lbrk = data.frame(let = let[1,lbrk], brk = lbrk, col = 
             founder.cols[match(let[1,lbrk], founder.cols[,1]),3],
             stringsAsFactors = F)
      rbrk = data.frame(let = let[2,rbrk], brk = rbrk, col =
             founder.cols[match(let[2,rbrk], founder.cols[,1]),3],
             stringsAsFactors = F)
      lbrk[,2] = snps[snps[,2] == names(chrlen)[c],3][lbrk[,2]]
      rbrk[,2] = snps[snps[,2] == names(chrlen)[c],3][rbrk[,2]]
      rect(xleft = rep(2 * c - offset, nrow(lbrk) + 1), ybottom = c(0, lbrk[,2]),
           xright = rep(2 * c, nrow(lbrk) + 1), ytop = c(lbrk[,2], chrlen[c]),
           col = lbrk[,3], border = -1)
      rect(xleft = rep(2 * c, nrow(rbrk) + 1), ybottom = c(0, rbrk[,2]),
           xright = rep(2 * c + offset, nrow(rbrk) + 1),
           ytop = c(rbrk[,2], chrlen[c]), col = rbrk[,3], border = -1)
    } # for(c)
  } # for(i)

  par(old.par)
} # plot.genotypes()

