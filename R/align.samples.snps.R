################################################################################
# Align all of the samples and SNPs in different variables. Remove samples and
# SNPs that do not occur in all variables and order the sample IDs and SNP IDs
# identically in all variables.
# This would be called before performing QTL mapping.
# Daniel Gatti
# Dan.Gatti@max.org
# Nov. 4, 2013
################################################################################
align.samples.snps = function(pheno, probs, K, addcovar, intcovar, snps) {
  if(is.null(rownames(pheno))) {
    stop("The samples IDs must be in rownames(pheno).")
  } # if(is.null(rownames(pheno))
  # Intersect the sample IDs for all variables with samples.
  samples = intersect(rownames(pheno), dimnames(probs)[[1]])
  if(!missing(K)) {
    if(is.null(rownames(K))) {
      stop("The samples IDs must be in rownames(K) and colnames(K).")
    } # if(is.null(rownames(K))
    samples = intersect(samples, rownames(K))
  } # if(!missing(K))
  if(!missing(addcovar)) {
    if(is.null(rownames(addcovar))) {
      stop("The samples IDs must be in rownames(addcovar).")
    } # if(is.null(rownames(addcovar))
    samples = intersect(samples, rownames(addcovar))  
  } # if(!missing(addcovar))
  if(!missing(intcovar)) {
    if(is.null(rownames(intcovar))) {
      stop("The samples IDs must be in rownames(intcovar).")
    } # if(is.null(rownames(intcovar))
    samples = intersect(samples, rownames(intcovar))
  } # if(!missing(intcovar))
  if(length(samples) == 0) {
    stop(paste("There were no intersecting sample IDs between all of",
         "the variables. Please verify that there are sample IDs in",
         "common in the dimnames of all of the variables."))
  } # if(length(samples) == 0)
  message(paste("Using", length(samples), "samples."))
  # Order all of the samples based on the order in pheno.
  pheno = pheno[rownames(pheno) %in% samples,,drop = FALSE]
  probs = probs[dimnames(probs)[[1]] %in% samples,,,drop = FALSE]
  probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,,drop = FALSE]
  if(!missing(K)) {
    K = K[rownames(K) %in% samples, colnames(K) %in% samples, drop = FALSE]
    K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), 
          colnames(K)), drop = FALSE]
  } # if(!missing(K))
  if(!missing(addcovar)) {
    addcovar = addcovar[rownames(addcovar) %in% samples,,drop = FALSE]
    addcovar = addcovar[match(rownames(pheno), rownames(addcovar)),,drop = FALSE]
  } # if(!missing(addcovar))
  if(!missing(intcovar)) {
    intcovar = intcovar[rownames(intcovar) %in% samples,,drop = FALSE]
    intcovar = intcovar[match(rownames(pheno), rownames(intcovar)),,drop = FALSE]
  } # if(!missing(intcovar))
  # Intersect the SNP IDs in SNPs and probs.
  snpids = intersect(snps[,1], dimnames(probs)[[3]])
  snps = snps[snps[,1] %in% snpids,,drop = FALSE]
  probs = probs[,,dimnames(probs)[[3]] %in% snpids,,drop = FALSE]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]]),drop = FALSE]
  if(length(snpids) == 0) {
    stop(paste("There were no intersecting SNP IDs between snps",
         "and probs. Please verify that there are SNP IDs in",
         "common in the snps[,1] and dimnames(probs)[[3]]."))
  } # if(length(samples) == 0)
  message(paste("Using", length(snpids), "markers."))
  return(list(pheno, probs, K, addcovar, intcovar, snps))
} # align.samples.snps()
