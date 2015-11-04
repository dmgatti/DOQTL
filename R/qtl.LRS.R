################################################################################
# Perform QTL mapping using the Likelihood Ratio Statistic (LRS).
# Allow a set of fixed covariates which will comprise the null model.
# The genotypes must be a multi-state probability array.
# This was originally created to perform mapping in the DO mice with 8
# founders.
# I calculate the LRS by taking the ratio of the residual variances
# multiplied by the number of samples.  This is equivalent to computing the 
# more complicated formulas.
# It produces the same values as:
#   ll1 = logLik(lm(y ~ fixed))
#   ll2 = logLik(lm(y ~ fixed + probs))
#   LRS = -2 * (ll1 - ll2)
# Tested this code with several phenotypes by comparing it to the code above.
# Daniel Gatti
# Dan.Gatti@jax.org
# July 1, 2011
################################################################################
# Arguments: pheno: numerical matrix with phenotypes to map.  These should be
#               transformed as needed (ie RankZ) and may not contain NA values.
#               Samples in rows and phenotypes in columns.  May also be a
#               vector with samples names in names.
#            probs: 3D numerical array with the genotype probabilities of the
#                  founders for each sample at each SNP. num.samples x 
#                  num.states x num.SNPs.  Each row of the matrix at each SNP
#                  (ie each sample) must sum to 1.
#            snps: data.frame with SNP IDs, chromosomes & bp locations in 
#                  columnds 1,2 & 3 respectively.
#            covar: binary matrix with 0/1 containing fixed covariates.
# Returns: list with 2 elements:
#          1. SNPs & LRS for the genotype at each SNP.
#          2. 3D array of coefficients for the full model at each SNP.
qtl.LRS = function(pheno, probs, snps, addcovar = NULL) {

  if(!is.matrix(pheno)) {
    pheno = as.matrix(pheno)
  } # if(!is.matrix(pheno))

  # Number of samples.
  n = nrow(pheno)

  # Number of SNPs.
  m = dim(probs)[3]

  # Number of phenotypes.
  ph = ncol(pheno)
  num.covar = 0
  null.var = NULL

  if(!is.null(addcovar)) {

    # Number of covariates.
    num.covar = ncol(addcovar)
    # Create the covariates matrix. 
    addcovar = as.matrix(cbind(Intercept = rep(1, n), addcovar, probs[,,1]))
    # Get the residual variance for the null model with just the fixed
    # covariates.
    null.var = 1.0 / diag(crossprod(qr.resid(qr(addcovar[,1:(num.covar + 1)]), pheno)))

  } else {

    # Create the covariates matrix. 
    addcovar = as.matrix(cbind(Intercept = rep(1, n), probs[,,1]))    
    # Get the residual variance for the null model with just the fixed
    # covariates.
    null.var = 1.0 / diag(crossprod(qr.resid(qr(addcovar[,1]), pheno)))

  } # else

  # Calculate the residual variances of the full model at each SNP.
  lrs = matrix(0, m, ph, dimnames = list(dimnames(probs)[[3]], colnames(pheno)))
  coef = array(0, c(m, ncol(addcovar), ph), dimnames = list(snps[,1], 
         colnames(addcovar), colnames(pheno)))

  # Find the SNPs with low MAF.
  maf = apply(apply(probs, 2, colSums), 1, min)
  run.pseudo = which(maf < 0.5)
  run.qr = which(maf >= 0.5)

  # First run the SNPs where we can use the QR decomposition.
  # The range of columns that contain the genotypes.
  rng = (ncol(addcovar) - dim(probs)[2] + 1):ncol(addcovar)
  for(s in run.qr) {
    addcovar[,rng] = probs[,,s]
    qrx = qr(addcovar)
    lrs[s,] = colSums(qr.resid(qrx, pheno)^2)
    coef[s,,] = qr.coef(qrx, pheno)
  } # for(s)

  # Then run the SNPs where one of the allele frequencies is too low and
  # we need to use the pseudo-inverse to solve the regression.
  for(s in run.pseudo) {
    addcovar[,rng] = probs[,,s]
    beta = pseudoinverse(t(addcovar) %*% addcovar) %*% t(addcovar) %*% pheno
    lrs[s,] = colSums((pheno - addcovar %*% beta)^2)
    coef[s,,] = beta
  } # for(s)

  lrs = -n * log(lrs * matrix(null.var, nrow(lrs), ncol(lrs), byrow = TRUE))
  p = pchisq(q = lrs, df = dim(probs)[[2]] - 1, lower.tail = FALSE)

  return(list(lrs = cbind(snps, lrs = lrs, lod = lrs / (2 * log(10)), p = p,
         neg.log10.p = -log(p, 10)), coef = coef))

} # qtl.LRS()


################################################################################
# Run permutations on the phenotypes to assess statistical significance of
# the QTL peaks.  Run all phenotypes together.  Be careful with this script.
# You must be able to permute the phenotypes and the fixed covariates 
# together for all phenotypes.
# Arguments: pheno: numerical matrix with phenotypes to map.  These should be
#               transformed as needed (ie RankZ) and may not contain NA values.
#               Samples in rows and phenotypes in columns.  May also be a
#               vector with sample names.
#            probs: 3D numerical array with the genotype probabilities of the
#                  founders for each sample at each SNP. num.samples x 
#                  num.states x num.SNPs.  Each row of the matrix at each SNP
#                  (ie each sample) must sum to 1.
#            snps: data.frame with SNP IDs, chromosomes & bp locations in 
#                  columnds 1,2 & 3 respectively.
#            covar: binary matrix with 0/1 containing fixed covariates.
#            nperm: integer, number of permutations to run.
#            return.val: character string containing either "lrs" or "p",
#                        indicating the type of return statistic.
permutations.qtl.LRS = function(pheno, probs, snps, addcovar, nperm = 1000,
                       return.val = c("lod", "p")) {

  return.val = match.arg(return.val)
  
  if(!is.matrix(pheno)) {
    pheno = as.matrix(pheno)
  } # if(!is.matrix(pheno))

  # Number of samples.
  n = nrow(pheno)

  # Number of SNPs.
  m = dim(probs)[3]

  # Number of phenotypes.
  ph = ncol(pheno)

  # Create the covariates matrix.
  null.var = NULL
  if(missing(addcovar)) {

    addcovar = as.matrix(cbind(Intercept = rep(1, n)))
    null.var = 1.0 / diag(crossprod(qr.resid(qr(addcovar[,1]), pheno)))
    addcovar = as.matrix(cbind(Intercept = rep(1, n), probs[,,1]))
    num.addcovar = 1

  } else {

    addcovar = as.matrix(addcovar)
    num.addcovar = ncol(addcovar)
    addcovar = as.matrix(cbind(Intercept = rep(1, n), addcovar))
    null.var = 1.0 / diag(crossprod(qr.resid(qr(addcovar[,1:(num.addcovar+1)]),
                          pheno)))
    addcovar = as.matrix(cbind(addcovar, probs[,,1]))

  } # else

  # Save the maximum LRS from each permutation.
  max.lrs = matrix(0, nperm, ph, dimnames = list(1:nperm, colnames(pheno)))
  lrs = matrix(0, m, ph, dimnames = list(dimnames(probs)[[3]], colnames(pheno)))

  # Get the sex of the samples.
  sex.col = grep("sex", colnames(addcovar), ignore.case = TRUE)
  sex = as.numeric(factor(addcovar[,sex.col])) - 1
  females = which(sex == 0)
  males   = which(sex == 1)

  # The range of columns that contain the genotypes.
  rng = (ncol(addcovar) - dim(probs)[2] + 1):ncol(addcovar)

  # Decide when to use QR or the pseudoinverse.
  rnk = rep(0, m)
  for(s in 1:m) {
    addcovar[,rng] = probs[,,s]
    rnk[s] = rank.condition(addcovar)$rank
  } # for(s)

  run.pseudo = which(rnk < ncol(addcovar))
  run.qr = which(rnk == ncol(addcovar))

  for(p in 1:nperm) {

    print(paste(p, "of", nperm))

    # Permute the phenotypes and fixed covariates.
    new.order = rep(0, nrow(pheno))
    new.order[females] = sample(females)
    new.order[males]   = sample(males)
    pheno = as.matrix(pheno[new.order,,drop = FALSE])
    addcovar = addcovar[new.order,,drop = FALSE]

    # First run the SNPs where we can use the QR decomposition.
    # The range of columns that contain the genotypes.
    for(s in run.qr) {
      addcovar[,rng] = probs[,,s]
      lrs[s,] = colSums(qr.resid(qr(addcovar), pheno)^2)
    } # for(s)

    # Then run the SNPs where one of the allele frequencies is too low and
    # we need to use the pseudo-inverse to solve the regression.
    for(s in run.pseudo) {
      addcovar[,rng] = probs[,,s]
      beta = pseudoinverse(t(addcovar) %*% addcovar) %*% t(addcovar) %*% pheno
      lrs[s,] = colSums((pheno - addcovar %*% beta)^2)
    } # for(s)
   
    # Note that we get the minimum residual variance to obtain the
    # minimum LRS.
    max.lrs[p,] = apply(lrs[,,drop = FALSE], 2, min)

  } # for(p)

  max.lrs = -n * log(max.lrs * matrix(null.var, nrow(max.lrs), 
             ncol(max.lrs), byrow = TRUE))

  if(return.val == "p") {

    max.lrs = pchisq(q = max.lrs, df = dim(probs)[[2]] - 1, lower.tail = FALSE)

  } # else

  return(max.lrs)

} # permutations.qtl.LRS()
