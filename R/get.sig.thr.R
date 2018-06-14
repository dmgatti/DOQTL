################################################################################
# Given a set of maximum LOD scores, calculate the significance thresholds
# for the given levels. Note that we use the method of Broman et.al., Genetics,
# 2006 to calculate the X chromosoem threshold.
# Daniel Gatti
# dan.gatti@jax.org
# Oct. 26, 2015
################################################################################
# Arguments: perms: either a vector with the maximum LOD score from permutations
#                   or, if Xchr == TRUE, a matrix with colnames 'A' and 'X' 
#                   containing the maximumum autosomal and X chr LOD scores from 
#                   permutations.
#            alpha: numeric vector, one or more numerical value that are the 
#                   alpha levels of the significane thresholds.
#            Xchr: logical indicating whether the X chromosome should have 
#                  a separate threshold. If TRUE, then perms must be a matrix
#                  with 2 named columns as described above.
# NOTE: If using p-values, don't forget to take -log10(p.value) of the 
#       permuations before calling this function.
get.sig.thr = function(perms, alpha = 0.05, Xchr = TRUE) {
  sig.thr = rep(0, length(alpha))
  if(Xchr) {
    if(!is.matrix(perms)) {
      stop(paste("'perms' is not a matrix. 'perms' must be a matrix",
           "with 2 columns, named 'A' and 'X'."))
    } # if(!is.matrix(perms))
    if(!(all(colnames(perms) %in% c("A", "X")))) {
      stop(paste("The colnames of 'perms' are not equal to 'A' and",
           "'X'. 'perms' must be a matrix, with 2 columns, named",
           "'A' and 'X'."))
    } # if(!(all(colnames(perms) %in% c("A", "X"))))
    chrlen = get.chr.lengths()
    len.auto = sum(chrlen[1:19])
    len.X = chrlen["X"]
    len.all = len.auto + len.X
    alpha.auto = 1.0 - (1.0 - alpha)^(len.auto / len.all)
    alpha.X    = 1.0 - (1.0 - alpha)^(len.X / len.all)
    sig.thr = cbind("A" = quantile(perms[,"A"], probs = 1.0 - alpha.auto, na.rm = TRUE),
                    "X" = quantile(perms[,"X"], probs = 1.0 - alpha.X, na.rm = TRUE))
    rownames(sig.thr) = alpha
  } else {
    sig.thr = quantile(perms, probs = 1.0 - alpha, na.rm = TRUE)
    names(sig.thr) = alpha
  } # else
  return(sig.thr)
} # get.sig.thr()
