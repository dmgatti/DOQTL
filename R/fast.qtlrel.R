################################################################################
# A faster impelmentation of the QTLRel algorithm for additive covariates with
# a kinship matrix.  For interactive covariates, use the main QTLRel function.
# This is only about twice as fast as QTLRel.
# Code adapted from Riyan Cheng's QTLRel package on CRAN.
# NOTE: No NAs are allowed anywhere! Filter your data before calling this.
# Daniel Gatti
# Dan.Gatti@jax.org
# Aug, 2, 2013
################################################################################
# Arguments: pheno: numeric vector with no NAs containing the phenotype.
#            probs:
fast.qtlrel = function(pheno, probs, K, addcovar, snps) {

  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = T) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) < 200)

  # Return value.
  retval = NULL

  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3], snp = snps[,1])
  class(prdat) = c(class(prdat), "addEff")
  err.cov = diag(length(pheno))
  if(!missing(K)) {
    K = as.matrix(K)
    vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
                EE = err.cov)
    vc = estVC(y = pheno, x = addcovar, v = vTmp)
    err.cov = matrix(0, nrow(K), ncol(K))
    num.cov = 1
    if(!missing(addcovar)) {
      num.cov = ncol(addcovar) + 1
    } # if(missing(addcovar))
    for(j in which(vc$nnl)) {
      err.cov = err.cov + vTmp[[j]] * vc$par[names(vTmp)[j]]
    } # for(j)
    rm(vTmp)
  } # if(!missing(K))

  # Invert the covariance matrix.
  eW = eigen(err.cov, symmetric = T)
  if (min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)) {
    stop("'W' is not positive definite")
  } else {
    eW$values[eW$values <= 0] = Inf
  } # else
  err.cov = diag(eW$values^-0.5) %*% t(eW$vector)
  rm(eW)

  # Add the intercept.
  if(missing(addcovar)) {
    addcovar = matrix(1, nrow = length(pheno), ncol = 1)
  } else {
    addcovar = as.matrix(cbind(rep(1, length(pheno)), addcovar))
  } # else
  colnames(addcovar)[1] = "Intercept"

  # Null model.
  ytmp = err.cov %*% pheno
  xtmp = err.cov %*% addcovar
  qr.null = qr(xtmp)
  ss.null = sum(qr.resid(qr.null, ytmp)^2)

  # Additive model for all SNPs.
  addx = cbind(addcovar, probs[,-1,1])
  ss = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), ncol(addx), dimnames = list(snps[,1],
         colnames(addx)))
  rng = (ncol(addcovar)+1):ncol(addx)
  for(s in 1:nrow(snps)) {
    addx[,rng] = probs[,-1,s]
    xtmp = err.cov %*% addx
    qr.add = qr(xtmp)
    ss[s] = sum(qr.resid(qr.add, ytmp)^2)
    coef[s,] = qr.coef(qr.add, ytmp)
  } # for(s)

  perc.var = 1.0 - (ss / ss.null)
  lrs = -length(pheno) * log(ss / ss.null)
  lod = lrs / (2 * log(10))

  return(list(lod = cbind(snps[,1:4], perc.var, lrs, lod),
         coef = coef))

} # fast.qtlrel()

