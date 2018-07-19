##########################################################################################
# Perform permutations of the linkage mapping model and return the maximum LOD (or minimum
# p-value) in each permutation on the autosomes and the X chromosome.
# Daniel Gatti
# dan.gatti@jax.org
# 
##########################################################################################
scanone.perm = function(pheno, pheno.col = 1, probs, addcovar, intcovar, snps,
               model = c("additive", "full"), path = ".", nperm = 1000,
               return.val = c("lod", "p")) {
  return.val = match.arg(return.val)
  if(!missing(intcovar)) {
    stop("Interactive covariates not yet implemented")
  } # if(!missing(intcovar))
  if(missing(addcovar)) {
    stop(paste("'addcovar' is required and must contain a variable called 'sex' in",
         "order to map correctly on the X chromosome. This is required even if all",
         "of your samples have the same sex."))
  } # if(missing(addcovar))
  if(length(grep("sex", colnames(addcovar), ignore.case = TRUE)) == 0) {
    stop(paste("'addcovar' is required and must contain a variable called 'sex' in",
         "order to map correctly on the X chromosome. This is required even if all",
         "of your samples have the same sex."))
  } # if(length(grep("sex",  ...
  if(is.null(rownames(addcovar))) {
    stop("rownames(addcovar) is null. The sample IDs must be in rownames(pheno).")
  } # if(is.null(rownames(addcovar)))
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
  snps  = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
  # Intersect common samples.
  samples = intersect(intersect(rownames(pheno), rownames(addcovar)), rownames(probs))
  pheno = pheno[samples,,drop = FALSE]
  addcovar = as.matrix(addcovar)
  addcovar = addcovar[samples,,drop = FALSE]
  probs = probs[samples,,,drop = FALSE]
  
  # Remove rows with NA, NaN or Inf.
  keep = which(rowMeans(is.na(addcovar) | is.nan(addcovar) | is.infinite(addcovar)) == 0 & 
               (is.na(pheno[,pheno.col]) |
               is.nan(pheno[,pheno.col]) | 
               is.infinite(pheno[,pheno.col])) == 0)
  if(length(keep) == 0) {
    stop(paste("There are no matching samples in pheno, addcovar and probs.",
        "Please verify that the sample IDs are in rownames(pheno),", 
        "rownames(addcovar) and rownames(probs) and that they all match."))
  } # if(any(dim(probs) == 0))
  # Match up the sample names in phenotype, addcovar and probs.
  pheno = pheno[keep,,drop = FALSE]
  probs = probs[keep,,,drop = FALSE]
  addcovar = addcovar[keep,, drop = FALSE]
  
  print(paste("Mapping with", nrow(pheno),"samples."))
  if(is.null(colnames(addcovar))) {
    colnames(addcovar) = paste("addcovar", 1:ncol(addcovar), sep = ".")
  } # if(is.null(colnames(addcovar)))
  # Get the autosomal and X chromosome markers to get separate thresholds on 
  # each chromosome.
  num.auto = get.num.auto(snps)
  auto.snps = which(snps[,2] %in% 1:num.auto)
  X.snps = which(snps[,2] == "X")
  # Create a separate set of probs for the females and males on Chr X.
  if(length(X.snps) > 0) {
    xprobs = array(0, c(nrow(probs), 2 * ncol(probs), length(X.snps)), 
             dimnames = list(rownames(probs), paste(rep(c("F", "M"), each = 8),
             LETTERS[1:8], sep = "."), snps[X.snps,1]))
    # Get the sex of the samples and place the probs in the correct columns.
    sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
    sex = as.numeric(factor(addcovar[,sex.col])) - 1
    females = which(sex == 0)
    males   = which(sex == 1)
    xprobs[females,1:8,] = probs[females,,X.snps]
    xprobs[males,9:16,]  = probs[males,,X.snps]
  } #   if(length(X.snps) > 0)
  # Perms is for one phenotype.
  perms = array(0, dim = c(nperm, length(pheno.col), 1), dimnames = 
          list(1:nperm, colnames(pheno)[pheno.col], 1))
  if(length(X.snps) > 0) {
    perms = array(0, c(nperm, length(pheno.col), 2), dimnames = 
            list(1:nperm, colnames(pheno)[pheno.col], c("A", "X")))
  } # if(length(X.snps) > 0)
  # Loop through each phenotype column.
  index = 1 # Index into the perms array.
  for(i in pheno.col) {
    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p) & !is.nan(p) & !is.infinite(p))
    auto.perms = permutations.qtl.LRS(pheno = p[keep], probs = probs[keep,,auto.snps], 
                 snps = snps[auto.snps,], addcovar = addcovar[keep,,drop = FALSE], 
                 nperm = nperm, return.val = return.val)
    X.perms = 0
    if(length(X.snps) > 0) {
      X.perms = permutations.qtl.LRS(pheno = p[keep], probs = xprobs[keep,,], 
                snps = snps[X.snps,], addcovar = addcovar[keep,,drop = FALSE], 
                nperm = nperm, return.val = return.val)
    } # if(length(X.snps) > 0)
    if(return.val == "lod" ) {
      auto.perms = auto.perms / (2 * log(10))
      if(length(X.snps) > 0) {
        X.perms = X.perms / (2 * log(10))
      } # if(length(X.snps) > 0)
    } # if(return.val == "lod" )
    if(length(X.snps) > 0) {
      perms[,index,] = cbind(auto.perms, X.perms)
    } else {
      perms[,index,1] = auto.perms
    } # if(length(X.snps) > 0)
    write.table(perms[,,index], file = paste(path, "/", colnames(pheno)[i], ".perms.txt",
          sep = ""), sep = "\t")
    index = index + 1
  } # for(i)
  return(perms)
  
} # scanone.perm()
