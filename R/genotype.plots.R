################################################################################
# Make a plot of DO genotypes on each chromosome.
# Feb. 14, 2011
# Daniel Gatti
# Dan.Gatti@jax.org
# Contains: get.chr.lengths, chr.skeletons, genomic.points, plot.genoprob,
#           write.genoprob.plots.
################################################################################
# Draw the skeleton of the chromosomes.
# Arguments: chr: vector with chr names to plot.
chr.skeletons = function(chr, chrlen = "mm10", ...) {
  
  if(chrlen[1] == "mm10") {
    chrlen = get.chr.lengths(genome = chrlen)
  } # if(chrlen[1] == "mm10")

  if(missing(chr)) {
    chr = names(chrlen)
  } # if(missing(chr))

  chrlen = chrlen[names(chrlen) %in% chr]

  # Plot a coordinate system to draw on.
  par(font = 2, font.axis = 2, font.lab = 2, las = 1, plt = c(0.05, 0.95, 
      0.1, 0.9))
  plot(-1, -1, col = 0, xlim = c(0, max(chrlen) + 10), 
       ylim = c(1, length(chrlen)),
       xaxs = "i", yaxt = "n", xlab = "", ylab = "", ...)

  # Draw chr skeletons
  for(i in 1:length(chrlen)) {
    lines(c(0, chrlen[i]), rep(length(chrlen) - i + 1, 2))
  } # for(i)

  # Draw light lines at 10 Mb and heavier lines at 50 Mb.
  abline(v = 0:30 * 10, col = "grey80")
  abline(v = 0:5  * 50, col = "grey50", lwd = 1.2)
  mtext(side = 2, line = 0.25, at = length(chrlen):1, text = names(chrlen))
  mtext(side = 1, line = 2.5, text = "Mb")

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
    points(chr, loc, pch = 3, cex = 2, lwd = 2, col = 2)

  } # if(!is.null(chr) & !is.null(loc))

} # genomic.points()


################################################################################
# Plot the genotypes using only the maximum probability at each SNP and arrange
# the colors to minimize jumping from one strand to the other.
# Arguments: x: numeric matrix, with smoothed probabilities. SNPs in rows, 
#                    genotypes in columns.
#            snps: matrix with SNP IDs in column 1, chromosomes in column 2 and
#                  Mb locations in column 3.
#            colors: data.frame with three columns containing the founder 
#                    letters, names and colors in columns 1, 2 & 3, respectively.
#                    Default = "DO", which automatically uses the CC founder 
#                    colors.
#            chrlen: named numeric vector containing chromosome lengths or "mm10"
#                    for the mouse chromosomes or "rn6" for rat.
#            type: Either "max", indicating that the two haplotypes with the
#                  maximum probability should be ploteed, or "probs", indicating
#                  that the actual haplotype probabilities should be plotted.
#            ...: other arguments to be passed to plot.
plot.genoprobs = function(x, snps, colors = "DO", chrlen = "mm10",
                          type = c("max", "probs"), legend = TRUE, ...) {

  if(any(rowSums(x) - 1 < -1e-8)) {
    stop(paste("plot.genoprobs: the probabilities in each row do not",
         "sum to 1."))
  } # if(!all(rowSums(x) - 1.0 < 1e-8))

  old.par = par(no.readonly = TRUE)
  states = colnames(x)

  if(max(snps[,3], na.rm = TRUE) > 300) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3], na.rm = TRUE) > 300)
  
  cross = attr(x, "cross")
  genome = chrlen
  type = match.arg(type)

  # Get the colors for the haplotype blocks and the genome to use.
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

    if(chrlen == "mm10") {
      chrlen = get.chr.lengths(genome = chrlen)
    } # if(missing(chrlen))

  } else {

    if(cross == "DO" | cross == "CC" | cross == "DOF1") {
      colors = do.colors
      chrlen = get.chr.lengths(genome = "mm10")
    } else if(cross == "HS") {
      colors = hs.colors
      chrlen = get.chr.lengths(genome = "mm10")
    } else if(cross == "HSrat") {
      colors = hsrat.colors
      chrlen = get.chr.lengths(genome = "rn6")
    } # else

  } # else

  # Subset the data to only include the SNPs in the snpfile.
  x = x[rownames(x) %in% snps[,1],]
  snps = snps[snps[,1] %in% rownames(x),]
  prmsth = x[match(snps[,1], rownames(x)),]

  # If we have 36 state diplotype probabilities, then condense them down to
  # 8 state haplotype probabilities.
  if(ncol(prsmth) > 8) {
    mat = get.diplotype2haplotype.matrix(colnames(prsmth))
    prsmth = prsmth %*% mat
    states = colnames(prsmth)
  } # if(ncol(prsmth) > 8)

  # Plot the chromosome skeletons.
  old.warn = options("warn")$warn
  options(warn = -1)
  char.chr = as.numeric(as.character(unique(snps[,2])))
  char.chr = which(is.na(char.chr))
  char.chr = as.character(unique(snps[,2]))[char.chr]
  chr = sort(as.numeric(as.character(unique(snps[,2]))), na.last = TRUE)
  options(warn = old.warn)
  chr[is.na(chr)] = char.chr

  par(plt = c(0.08, 0.95, 0.05, 0.88))
  chr.skeletons(chr, chrlen = chrlen, ...)
  chrlen = chrlen[names(chrlen) %in% snps[,2]]

  # If we have 36 state probs, convert them to 8 state.
  if(ncol(x) == 36) {

    mat = get.diplotype2haplotype.matrix(colnames(x))
    x = x %*% mat

  } # if(ncol(x) == 36)

  if(type == "max") {
    plot.genoprobs.max(x = x, snps = snps, colors = colors, chrlen = chrlen,
                       offset = 0.35)
  } else if(type == "probs") {
    plot.genoprobs.probs(x = x, snps = snps, colors = colors, chrlen = chrlen,
                         offset = 0.6)
  } else {
    stop(paste("plot.genoprobs: Invalid plot type", type))
  } # else

  usr = par("usr")
  if(legend) {
    legend(x = 160, y = 10, legend = colors[,2], col = colors[,3], pch = 15,
           bg = "white")
  } # if(legend)

} # plot.genoprobs()


plot.genoprobs.probs = function(x, snps, colors = "DO", chrlen = "mm10",
                                offset) {

  # Synch up the SNPs and probs.
  x = x[snps[,1],]
  stopifnot(all(rownames(x) == snps[,1]))

  # Make a matrix with the start and end of each prob.
  pltmat = cbind(0, x)
  pltmat = (apply(pltmat, 1, cumsum) - 0.5) * offset

  # Split the data by chormosome.
  pltmat = split(t(pltmat), snps[,2])
  # Order by chromosome.
  pltmat = pltmat[order(as.numeric(names(pltmat)))]
  # Remove chromsomes without data.
  pltmat = pltmat[sapply(pltmat, length) > 0]
  # Reform the data into matrices.
  pltmat = lapply(pltmat, function(z) { matrix(z, ncol = 9) })
  pltmat = lapply(pltmat, t)

  pos = split(snps[,3], snps[,2])
  pos = pos[names(pltmat)]

  for(chr in 1:length(pltmat)) {

    chrmid = length(chrlen) - chr + 1

    for(i in 1:ncol(pltmat[[chr]])) {

      rect(xleft = pos[[chr]][i], ybottom = pltmat[[chr]][1:8,i] + chrmid, 
           xright = pos[[chr]][i+1], ytop = pltmat[[chr]][2:9,i] + chrmid,
           col = colors[,3], border = colors[,3])

     } # for(i)

  } # for(chr)

} # plot.genoprobs.probs()


plot.genoprobs.max = function(x, snps, colors, chrlen, offset) {

  if(!is.matrix(x)) {
    stop("call.haps: x must be a matrix.")
  } # if(!is.matrix(x))

  if(ncol(x) > 8) {
    stop("call.haps: x must have only 8 columns.")
  } # if(ncol(x) > 8)

  # Get the unique chromosomes.
  chr = unique(snps[,2])
  chr = as.character(chr)
  old.warn = options("warn")$warn
  options(warn = -1)
  na.chr = is.na(as.numeric(chr))
  chr = factor(chr, levels = c(sort(as.numeric(chr[!na.chr])), chr[na.chr]))
  options(warn = old.warn)

  # Split the data by chromosome.
  chr = factor(snps[,2], levels = chr)
  snps = split(snps, chr)
  x    = split(data.frame(x), chr)
  stopifnot(length(snps) == length(x))

  # Loop through each chromosome.
  for(i in 1:length(x)) {

    # Get the breakpoints.
    # Heterozygous breakpoints. Round 2 * probs.
    x[[i]] = as.matrix(x[[i]])
    pr = round(x[[i]] * 2)

    # Check for rows that do not sum to 2 (for 2 strands).
    wh = which(rowSums(pr) != 2)
    if(length(wh) > 0) {
      for(j in wh) {
        tmp = sort(x[[i]][j,])
        pr[j,] = colnames(pr) %in% names(tmp)[7:8]
      } # for(wh)
    } # if(length(wh) > 0)

    pr[nrow(pr),] = 0
    brk = diff(rbind(0, pr > 0))
    st  = apply(brk > 0, 2, which)
    en  = apply(brk < 0, 2, which)

    # Apply doesn't have 'simplify' argument, so we have to check for this.
    if(is.matrix(st)) {

      st = data.frame(st)
      en = data.frame(en)

    } # if(is.matrix(st))

    # Homozygous breakpoints.
    brk2 = diff(rbind(0, pr > 1))
    if(max(brk2) > 0) {
      st = mapply(c, st, apply(brk2 > 0, 2, which), SIMPLIFY = FALSE)
      en = mapply(c, en, apply(brk2 < 0, 2, which), SIMPLIFY = FALSE)
    } # sum(brk2)

    # Combine start and end into one matrix.
    stopifnot(sapply(st, length) == sapply(en, length))
    st.en = cbind(unlist(st), unlist(en))

    # This puts the founder IDs in rownames.
    rownames(st.en) = substring(rownames(st.en), 1, 1)
    st.en = st.en[order(st.en[,1]),]

    # Get the left and right boundaries on the (arbitrary) plus strand.
    num.breaks = 500  # Might be too small for one chr in advanced intercrosses.
    left  = rep(0, num.breaks)
    right = rep(0, num.breaks)
    col   = rep("", num.breaks)
    j = 1
    idx = 1
    while(idx < nrow(st.en)) { 

      left[j]   = snps[[i]][st.en[idx, 1], 3]
      next.step = which(st.en[,1] == st.en[idx,2])[1]
      if(is.na(next.step)) {
        right[j] = snps[[i]][st.en[idx, 2], 3]
        col[j]   = colors[which(colors[,1] == rownames(st.en)[idx]),3]
        st.en = st.en[-idx,,drop = FALSE]
        break;
      } else {
        right[j] = snps[[i]][st.en[next.step, 2], 3]
        col[j]   = colors[which(colors[,1] == rownames(st.en)[idx]),3]
        st.en = st.en[-idx,,drop = FALSE]
        idx = next.step - 1
        j = j + 1
      } # else

    } # while(idx < nrow(st.en))

    remove = which(left == 0)
    left  = left[-remove]
    right = right[-remove]
    col   = col[-remove]

    # Draw the plus strand rectangles.
    rect(left, length(chrlen) - i + offset + 1, right, length(chrlen) - i + 1,
         col = col, border = NA)

    if(length(st.en) > 0) {

      # Get the left and right boundaries on the (arbitrary) minus strand.
      left  = snps[[i]][st.en[,1],3]
      right = snps[[i]][st.en[,2],3]
      col   = colors[match(rownames(st.en), colors[,1]), 3]

      # Draw the minus strand rectangles.
      rect(left, length(chrlen) - i - offset + 1, right, length(chrlen) - i + 1,
           col = col, border = NA)

    } # if(nrow(st.en) > 0)

  } # for(i)

} # plot.genoprobs.max()

################################################################################
# Make genotype plots for all of the files in a directory.
# This requires that the *.Rdata files have been made.
write.genoprob.plots = function(path = ".", snps, type = c("max", "probs")) {

  type = match.arg(type)
  files = dir(path, pattern = "genotype.probs.Rdata", full.names = TRUE)
  
  if(!is.null(files)) {
  
    prsmth = NULL
    for(f in files) {
    
      print(f)
      
      # This loads in 'prsmth'.
      load(f)
      sample = gsub(paste("^", path, "/|\\.genotype\\.probs\\.Rdata$", sep = ""),
                    "", f)
                    
      png(paste(sample, "genotype.probs.png", sep = "."), width = 1000,
          height = 1000, res = 144)
      plot.genoprobs(x = prsmth, snps = snps, type = type, main = sample)
      dev.off()
      
    } # for(f)
    
  } # if(!is.null(files))
  
} # write.genoprob.plots()
