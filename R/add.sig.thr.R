################################################################################
# Given a set of significance thresholds and chr lengths that have been used
# in a plot, draw the significance lines.
# NOTE: A QTL PLOT MUST HAVE BEEN DRAWN ALREADY.
# Daniel Gatti
# dan.gatti@jax.org
# Oct. 26, 2015
################################################################################
add.sig.thr = function(sig.thr, sig.col = "red", chrsum) {

  # If sig.thr is a matrix, then we have separate autosome and X thresholds.
  if(is.matrix(sig.thr)) {

    # Make sure that the thresholds and colors have the same length.
    if(length(sig.col) < nrow(sig.thr)) {
      stop(paste("The number of colors is less than the number of thresholds.",
           "Please provide as many colors as there are rows in the sig.thr matrix."))
    } # if(length(sig.col) < nrow(sig.thr))

    # Get the coordinates to the end of the autosomes.
    old.warn = options("warn")$warn
    options(warn = -1)    
    tmp = as.numeric(names(chrsum))
    auto.len = max(chrsum[!is.na(tmp)])
    X.len = max(chrsum[is.na(tmp)])
    options(warn = old.warn)

    for(i in 1:nrow(sig.thr)) {
      lines(x = c(0, auto.len), y = rep(sig.thr[i,"A"], 2), lwd = 2,
            col = sig.col[i])
      lines(x = c(auto.len, X.len), y = rep(sig.thr[i,"X"], 2), lwd = 2,
             col = sig.col[i])
    } # for(i)

  } else {

    # Make sure that the thresholds and colors have the same length.
    if(length(sig.col) < length(sig.thr)) {
      stop(paste("The number of colors is less than the number of thresholds.",
           "Please provide as many colors as there are thresholds in sig.thr."))
    } # if(length(sig.col) < nrow(sig.thr))

    abline(h = sig.thr, lwd = 2, col = sig.col)

  } # else

} # add.sig.thr
