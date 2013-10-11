################################################################################
# Given a set of mapping statistics and a set of SNPs, make a simple QTL plot.
# Arguments: doqtl: list containing two elements from scanone.
#                   lod: data.frame with chr locations and LOD score.
#                   coef: matrix with model coeffieints.
#            stat.name: character, the name of the mapping statistic (used as
#                       y-axis label.
#            sig.thr: vector, numeric, a set of significane thresholds to plot.
#            sig.col: vector, color, a set of colors to use for each
#                     significance threshold. Must be same length as sig.thr.
# To plot a subset of chromosomes, feed in a subset of SNPs.
# Get the chr lengths from the org.Mm.eg.db pacakge.
qtl.plot = function(doqtl, stat.name = "LOD", sig.thr = NULL, sig.col = "red", ...) {
  
  chrlen = get.chr.lengths()

  # Get the call and arguments.
  call = match.call()

  lod = doqtl$lod$A
  if(any(names(doqtl$lod) == "X")) {
    lod = rbind(doqtl$lod$A, doqtl$lod$X)
  } # if(length(lod) > 1)

  # Get chr lengths and locations.  Create an x-axis based on Genome Mb.
  if(max(lod[,3], na.rm = T) > 200) {
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
  chrlen = cumsum(chrlen)
  chrlen = c(0, chrlen)
  chrmid = chrlen[-length(chrlen)] + (diff(chrlen) * 0.5)

  # Add the preceeding chromosome lengths to each SNP position.
  for(c in 2:length(unique.chr)) {
    rows = which(lod[,2] == unique.chr[c])
    gmb[rows] = gmb[rows] + chrlen[c]
  } # for(c)

  # Make the basic plot.
  par(font = 2, font.lab = 2, font.axis = 2, las = 1, xaxs = "i")
  if("ylim" %in% names(call)) {
    plot(gmb, lod[,7], col = 0, xlab = "Genome Mb", ylab = stat.name, ...)
  } else {
    if(!missing(sig.thr)) {
    plot(gmb, lod[,7], col = 0, xlab = "Genome Mb", ylab = stat.name,
         ylim = c(0, max(lod[,7], sig.thr, na.rm = T) * 1.05), ...)
    } else {
      plot(gmb, lod[,7], col = 0, xlab = "Genome Mb", ylab = stat.name,
           ylim = c(0, max(lod[,7], na.rm = T) * 1.05), ...)
    } # else
  } # else

  lod = cbind(lod, gmb)
  lod = split(lod, lod[,2])

  usr = par("usr")
  rect(chrlen[2 * 1:(length(chrlen) * 0.5) - 1], usr[3],
       chrlen[2 * 1:(length(chrlen) * 0.5)], usr[4], col = rgb(0.8,0.8,0.8),
       border = NA)
  rect(usr[1], usr[3], usr[2], usr[4], border = 1)
  text(chrmid, 0.95 * usr[4], names(chrlen)[-1])
  lapply(lod, function(z) { points(z[,8], z[,7], type = "l", lwd = 2)})
  if(!is.null(sig.thr)) {
    abline(h = sig.thr, col = sig.col)
  } # if(!is.null(sig.thr))

} # qtl.plot()

