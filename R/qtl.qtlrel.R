################################################################################
# QTLRel multi-founder linear model script.
# Daniel Gatti
# Dan.Gatti@jax.org
# June 15, 2011
# Jan. 11, 1012: Added code to extract the proportion of variance explained by
#                each QTL.  Riyan sent me a new QTLRel version (0.2.11).
# Aug. 1, 2013: Separated X Chr mapping and added interactive covariates.
################################################################################
# Arguments: pheno: numeric vector with phenotype values, sample IDs in names.
#                   No NAs are allowed.
#            probs: numeric 3D array with samples in dim[1], states in dim[2],
#                  and SNPs in dim[3]. No NAs are allowed.
#            K: numeric matrix that is the QTLRel additive genetic matrix 
#               produces using SNPs and genMatrix(). No NAs are allowed.
#            addcovar: numeric matrix with additive covariates (excluding the 
#                      intercept). No NAs are allowed.
#            intcovar: numeric matrix with interactive covariates (excluding the 
#                      intercept). No NAs are allowed.
#            snps: snps to use in the analysis. nrow must equal dim(probs)[3].
# Returns: LOD from QTL mapping with SNPs.
qtl.qtlrel = function(pheno, probs, K, addcovar, intcovar, snps) {
  pheno = as.matrix(pheno)
  if(!missing(K)) {
    K = as.matrix(K)
  } # if(!missing(K))
  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) > 200)

  # Return value.
  retval = NULL

  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3],
               snp = snps[,1])
  vTmp = list(AA = NULL, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
              EE = diag(length(pheno)))

  if(!missing(K)) {
    vTmp$AA = 2 * K
  } # if(!missing(K))

  # This tells QTLRel to fit the additive model.
  class(prdat) = c(class(prdat), "addEff") 
  vc = NULL
  if(missing(addcovar)) {
    # No covariates.
    vc = estVC(y = pheno, v = vTmp)
    res = scanOne(y = pheno, prdat = prdat, vc = vc, test = "None", numGeno = TRUE)
  } else {
    vc = estVC(y = pheno, x = addcovar, v = vTmp)
    if(missing(intcovar)) {
      # Additive covariates only.
      res = scanOne(y = pheno, x = addcovar, prdat = prdat, vc = vc, 
            numGeno = TRUE, test = "None")
    } else {
      # Additive and interactive covariates.
      res = scanOne(y = pheno, x = addcovar, prdat = prdat, vc = vc, 
            intcovar = intcovar, numGeno = TRUE, test = "None")
    } # else
  } # else

  # Convert the model coefficients to a matrix.
  coef = matrix(unlist(res$parameters), length(res$parameters),
         length(res$parameters[[1]]), dimnames = list(res$snp,
         names(res$parameters[[1]])), byrow = TRUE)

  # Return the LRS, LOD, p-value and -log10(p-value).
  p = pchisq(q = res$p, df = dim(probs)[[2]] - 1, lower.tail = FALSE)
  return(list(lod = cbind(snps[,1:4], perc.var = res$v, lrs = res$p,
           lod = res$p / (2 * log(10)), p = p,
           neg.log10.p = -log(p, 10)), coef = coef))

} # qtl.qtlrel()
