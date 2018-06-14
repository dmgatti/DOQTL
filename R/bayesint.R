####################################################################################
# Find a Bayesian Credible Interval for a QTL, given the LOD score and Chr
# locations.  We do this by scaling the area under the LOD score to 1 and finding
# 95% of the area, centered on the maximum QTL.
# Daniel Gatti
# Dan.Gatti@jax.org
# May 15, 2012
####################################################################################
# Arguments: qtl: QTL object as produced by scanone.
#            chr: character, chromosome ID.
#            prob: numeric, value between 0 and 1.
#            expandtomarkers: boolean, expand the QTL interval to the outermost
#                 flanking markers if TRUE.  Default = FALSE.
bayesint = function(qtl, chr, prob = 0.95, expandtomarkers = TRUE) {
  if(missing(qtl)) {
    stop("bayesint: The qtl cannot be null. Please supply a QTL object.")
  } # if(is.null(qtl))
  
  if(missing(chr)) {
    stop(paste("bayesint: The chromosome cannot be null."))
  } else if(!chr %in% c(1:19, "X")) {
    stop(paste("bayesint: The chromosome must be 1 to 19 or X."))
  } # else if
  
  if(prob < 0 | prob > 1) {
   stop(paste("bayesint: The probability must between 0 and 1."))
  } # if(prob < 0 | prob > 1)
  
  old.warn = options("warn")$warn
  options(warn = -1)
  if(!is.na(as.numeric(chr))) {
    qtl = qtl$lod$A
  } else {
    qtl = qtl$lod$X
  } # else
  options(warn = old.warn)
  
  qtl[,1] = as.character(qtl[,1])
  qtl[,2] = as.character(qtl[,2])
  qtl[,3] = as.numeric(qtl[,3])
  qtl[,7] = as.numeric(qtl[,7])
  qtl = qtl[qtl[,2] == chr,]
  pos = qtl[,3]
  
  if(any(is.na(pos))) {
    remove = which(is.na(pos))
    qtl = qtl[-remove,]
    pos = pos[-remove]
  } # if(any(is.na(pos)))
  
  # Make a set of 10,000 intervals so that we can integrate numerically.
  breaks = approx(x = pos, y = 10^qtl[,7], xout = seq(pos[1], pos[length(pos)],
                  length.out = 1e5))
  widths  = diff(breaks$x)
  heights = breaks$y[-1] + breaks$y[-length(breaks$y)]
  trapezoids = 0.5 * heights * widths
  # Normalize the area to 1.
  trapezoids = trapezoids / sum(trapezoids)
  # This code copied from the R/qtl bayesint() function by Karl Broman.
  ord = order(breaks$y[-length(breaks$y)], decreasing = TRUE)
  wh  = min(which(cumsum(trapezoids[ord]) >= prob))
  int = range(ord[1:wh])
  # Find the left & right SNP.
  left.snp  = c(NA, qtl[1,2], breaks$x[int][1], 
                approx(qtl[,3], qtl[,4], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,5], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,6], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,7], breaks$x[int][1])$y)
  max.snp   = qtl[which.max(qtl[,7]),]
  right.snp = c(NA, qtl[1,2], breaks$x[int][2], 
                approx(qtl[,3], qtl[,4], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,5], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,6], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,7], breaks$x[int][2])$y)
  if(expandtomarkers) {
    # Find the left & right SNP.
    left.snp  = qtl[max(which(breaks$x[int][1] >= qtl[,3])),]
    max.snp   = qtl[which.max(qtl[,7]),]
    right.snp = qtl[min(which(breaks$x[int][2] <= qtl[,3])),] 
  } # if(expandtomarkers)
  retval = rbind(left.snp, max.snp, right.snp)
  retval[,3] = round(as.numeric(retval[,3]), digits = 6)
  retval[,4] = round(as.numeric(retval[,4]), digits = 6)
  retval[,5] = round(as.numeric(retval[,5]), digits = 6)
  retval[,6] = round(as.numeric(retval[,6]), digits = 6)
  retval$lod = as.numeric(retval[,7])
  return(retval)
} # bayesint()
