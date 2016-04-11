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
#            probs: 3D numeric array of haplotype probabilities.
#            K: numeric matrix.
fast.qtlrel = function(pheno, probs, K, addcovar, snps) {

  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) < 200)

  err.cov = NULL

  if(!missing(K)) {

    K = as.matrix(K)
    mod = NULL
    # Force the variance component estimates to be positive.
    if(missing(addcovar)) {
      mod = regress(pheno ~ 1, ~K, pos = c(TRUE, TRUE))
    } else {
      mod = regress(pheno ~ addcovar, ~K, pos = c(TRUE, TRUE))
    } # else

    err.cov = mod$sigma[1] * K + mod$sigma[2] * diag(length(pheno))

    # Invert the covariance matrix.
    eW = eigen(err.cov, symmetric = TRUE)
    if (min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)) {
      stop("fast.qtlrel: W is not positive definite")
    } else {
      eW$values[eW$values <= 0] = Inf
    } # else
    err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
    rm(eW)

  } # if(!missing(K))

  # Add the intercept.
  if(missing(addcovar)) {
    addcovar = matrix(1, nrow = length(pheno), ncol = 1, dimnames =
               list(rownames(pheno), "Intercept"))
  } else {
    addcovar = as.matrix(cbind(rep(1, length(pheno)), addcovar))
    colnames(addcovar)[1] = "Intercept"
  } # else

  # Remove strain A as the basis.
  probs = probs[,-1,]

  # Null model.
  ss.null = 0
  if(!is.null(err.cov)) {

    ytmp = err.cov %*% pheno
    xtmp = err.cov %*% addcovar
    qr.null = qr(xtmp)
    ss.null = sum(qr.resid(qr.null, ytmp)^2)

  } else {

    qr.null = qr(addcovar)
    ss.null = sum(qr.resid(qr.null, pheno)^2)

  } # else

  # Additive model for all SNPs.
  addx = cbind(addcovar, probs[,,1])
  ss = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), ncol(addx), dimnames = list(snps[,1],
         colnames(addx)))
  rng = (ncol(addcovar)+1):ncol(addx)

  perc.var = 0
  lrs = 0
  lod = 0

  # Find the SNPs with low MAF.
  maf = apply(apply(probs, 2, colSums), 1, min)
  run.pseudo = which(maf < 0.5)
  run.qr = which(maf >= 0.5)

  if(!is.null(err.cov)) {

    # First run the SNPs where we can use the QR decomposition.
    for(s in run.qr) {
      addx[,rng] = probs[,,s]
      xtmp = err.cov %*% addx
      qr.add = qr(xtmp)
      ss[s] = sum(qr.resid(qr.add, ytmp)^2)
      coef[s,] = qr.coef(qr.add, ytmp)
    } # for(s)

    # Then run the SNPs where one of the allele frequencies is too low and
    # we need to use the pseudo-inverse to solve the regression.
    for(s in run.pseudo) {
      addx[,rng] = probs[,,s]
      xtmp = err.cov %*% addx
      beta = pseudoinverse(t(xtmp) %*% xtmp) %*% t(xtmp) %*% ytmp
      ss[s] = colSums((ytmp - xtmp %*% beta)^2)
      coef[s,] = beta
    } # for(s)

  } else {

    # First run the SNPs where we can use the QR decomposition.
    for(s in run.qr) {

      addx[,rng] = probs[,,s]
      qr.add = qr(addx)
      ss[s] = sum(qr.resid(qr.add, pheno)^2)
      coef[s,] = qr.coef(qr.add, pheno)

    } # for(s)

    # Then run the SNPs where one of the allele frequencies is too low and
    # we need to use the pseudo-inverse to solve the regression.
    for(s in run.pseudo) {
      addx[,rng] = probs[,,s]
      beta = pseudoinverse(t(addx) %*% addx) %*% t(addx) %*% pheno
      ss[s] = colSums((pheno - addx %*% beta)^2)
      coef[s,] = beta
    } # for(s)

  } # else

  perc.var = 100 * (1.0 - (ss / ss.null))
  lrs = -length(pheno) * log(ss / ss.null)
  lod = lrs / (2 * log(10))

  # Get the p-value from the LRS.
  p = pchisq(q = -length(pheno) * log(ss / ss.null), 
      df = ncol(addx) - ncol(addcovar), lower.tail = FALSE)

  return(list(lod = cbind(snps[,1:4], perc.var  = perc.var, lrs = lrs, 
         lod = lod, p = p, neg.log10.p = -log(p, 10)), coef = coef))

} # fast.qtlrel()

