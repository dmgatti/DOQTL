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
#            test: character, One of "None", "F", or "Chisq" passed down to
#                  the QTLRel scanOne function.  "None" = LOD.
# Returns: LOD from QTL mapping with SNPs.
qtl.qtlrel = function(pheno, probs, K, addcovar, intcovar, snps, test = "None") {

  pheno = as.matrix(pheno)
  K = as.matrix(K)

  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = T) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) < 200)

  # Return value.
  retval = NULL

  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3],
               snp = snps[,1])
  vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
              EE = diag(length(pheno)))

  # This tells QTLRel to fit the additive model.
  class(prdat) = c(class(prdat), "addEff") 
  vc = NULL
  if(missing(addcovar)) {
    # No covariates.
    vc = estVC(y = pheno, v = vTmp)
    res = scanOne(y = pheno, prdat = prdat, vc = vc, test = test, numGeno = T)
  } else {
    vc = estVC(y = pheno, x = addcovar, v = vTmp)
    if(missing(intcovar)) {
      # Additive covariates only.
      res = scanOne(y = pheno, x = addcovar, prdat = prdat, vc = vc, 
            numGeno = T, test = test)
    } else {
      # Additive and interactive covariates.
      res = scanOne(y = pheno, x = addcovar, prdat = prdat, vc = vc, 
            intcovar = intcovar, numGeno = T, test = test)
    } # else
  } # else

  # Convert the model coefficients to a matrix.
  coef = matrix(unlist(res$parameters), length(res$parameters),
         length(res$parameters[[1]]), dimnames = list(res$snp,
         names(res$parameters[[1]])), byrow = T)

  return(list(lod = cbind(snps[,1:4], perc.var = res$v, lrs = res$p,
         lod = res$p / (2 * log(10))), coef = coef))

} # qtl.qtlrel()

