scanone.perm = function(pheno, pheno.col = 1, probs, addcovar, intcovar, snps,
               model = c("additive", "full"), path = ".", nperm = 1000) {

  if(!missing(intcovar)) {
    stop("Interactive covariates not yet implemented")
  } # if(!missing(intcovar))

  if(is.null(rownames(pheno))) {
    stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
  } # if(is.null(rownames(pheno)))

  if(is.character(pheno.col)) {
    pheno.col = match(pheno.col, colnames(pheno))
  } # if(is.character(pheno.col))
  
  if(!file.exists(path)) {
    stop(paste("The path", path, "does not exist."))
  } # if(!file.exists(path))

  probs = filter.geno.probs(probs)

  # Match up the sample names in the phenotype and probstype matrices.
  pheno = pheno[rownames(pheno) %in% dimnames(probs)[[1]],,drop = FALSE]
  probs = probs[dimnames(probs)[[1]] %in% rownames(pheno),,]
  probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,]
  probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]

  if(any(dim(probs) == 0)) {
    stop(paste("There are no matching samples in pheno and probs. Please",
         "verify that the sample IDs in rownames(pheno) match the sample",
         "IDs in dimnames(probs)[[1]]."))
  } # if(any(dim(probs) == 0))

  probs = filter.geno.probs(probs)
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
  
  if(!missing(addcovar)) {
    addcovar = as.matrix(addcovar)
    addcovar = addcovar[rowMeans(is.na(addcovar)) == 0,,drop = FALSE]
    pheno = pheno[rownames(pheno) %in% rownames(addcovar),]
    addcovar = addcovar[rownames(addcovar) %in% rownames(pheno),,drop = FALSE]
    addcovar = addcovar[match(rownames(pheno), rownames(addcovar)),,drop = FALSE]
    probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]

    if(is.null(colnames(addcovar))) {
      colnames(addcovar) = paste("addcovar", 1:ncol(addcovar), sep = ".")
    } # if(is.null(colnames(addcovar)))
  } # if(!is.null(addcovar))

  if(!missing(intcovar)) {
    covar = as.matrix(intcovar)
    intcovar = intcovar[rowMeans(is.na(intcovar)) == 0,,drop = FALSE]
    pheno = pheno[rownames(pheno) %in% rownames(intcovar),]
    intcovar = intcovar[rownames(intcovar) %in% rownames(pheno),,drop = FALSE]
    intcovar = intcovar[match(rownames(pheno), rownames(intcovar)),,drop = FALSE]
    probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]

    if(is.null(colnames(intcovar))) {
      colnames(intcovar) = paste("intcovar", 1:ncol(intcovar), sep = ".")
    } # if(is.null(colnames(intcovar)))
  } # if(!is.null(intcovar))

  for(i in pheno.col) {
    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p))
    if(!missing(addcovar)) {
      perms = permutations.qtl.LRS(pheno = p[keep], probs = probs[keep,,], 
              snps = snps, addcovar = addcovar[keep,,drop=FALSE], nperm = nperm)
    } else {
      perms = permutations.qtl.LRS(pheno = p[keep], probs = probs[keep,,], 
              snps = snps, nperm = nperm)
    } # else
    perms = sort(perms / (2 * log(10)))
    write(perms, paste(path, "/", colnames(pheno)[i], ".perms.txt", sep = ""),
          sep = "\t")
  } # for(i)
  return(perms)
  
} # scanone.perm()
