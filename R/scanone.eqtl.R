################################################################################
# Use the methods in Matrix eQTL to perform fast eQTL mapping for expression
# data. Method is from:
# Shabalin AA., Matrix eQTL: ultra fast eQTL analysis via large matrix
# operations, Bioinformatics. 2012 May 15;28(10):1353-8. PMID: 22492648
# Andrey Shabalin collaborated and wrote the code.
# Daniel Gatti
# Dan.Gatti@jax.org
# May 22, 2012
################################################################################
# Arguments: expr: numeric matrix with samples in rows and genes in columns.
#                  Rows and columns must contain sample and gene names.
#            probs: numeric array with samples in dim[1], 8 founders in dim[2],
#                   SNPs in dim[3]. Rows, columns and slices must have names.
#            K: numeric matrix with kinship coefficients. Samples in rows
#                 and columns. Rows and columns must have names.
#            addcovar: numeric matrix with additive covariates. samples in 
#                   rows, covariates in columns. Do not include the intercept.
#                   Rows and columns must have names.
# Returns: numeric matrix with LOD values for each gene and snp.
scanone.eqtl = function(expr, probs, K, addcovar, snps, sex) {
  time1 = proc.time()[3]
  # expr, probs and snps cannot be null.
  if(missing(expr) || missing(probs)) {
    stop(paste("One of the required arguments is null.  Please make sure",
         "that expr and probs are not null."))
  } # if(is.null(expr) | ...
  expr = t(expr)
  # Verify that the number of rows and columns match in the data.
  # Same number of samples in expr, probs and addcovar.
  if(ncol(expr) != dim(probs)[[1]]) {
    stop(paste("The number of columns in the expression matrix (",
         ncol(expr),") is not equal to the number of rows in probs (",
         dim(probs)[[1]],")."))
  } # if(ncol(expr) != dim(probs)[[1]])
  
  # Verify that the number of samples is the same in addcovar and expr.
  if(!missing(addcovar)) {
    if(nrow(addcovar) != ncol(expr)) {
      stop(paste("The number of rows in addcovar (", nrow(addcovar),") must",
           "equal the number of columns in expr (", nrow(expr) ,")."))
    } # if(nrow(addcovar) != nrow(expr))
  } # if(!missing(addcovar))
  # Verify that the number of samples is the same in the kinshp and 
  # expression matrices.
  if(!missing(K)) {
    if(nrow(K) != ncol(expr)) {
      paste(stop("The number of rows in the kinship matrix (", nrow(K),
            ") does not equal the number of columns in expr (", ncol(expr),
            ")."))
    } # if(nrow(K) != ncol(expr))
    if(ncol(K) != ncol(expr)) {
      paste(stop("The number of columns in the kinship matrix (", ncol(K),
            ") does not equal the number of columns in expr (", ncol(expr),
            ")."))
    } # if(ncol(K) != ncol(expr))
  } # if(!missing(K))
  # Match up the sample IDs.
  if(is.null(colnames(expr))) {
    stop(paste("expr must have sample IDs in the column names."))
  } # if(is.null(colnames(expr)))
  probs = probs[match(dimnames(probs)[[1]], colnames(expr)),,]
  if(!missing(addcovar)) {
    addcovar = addcovar[match(rownames(addcovar), colnames(expr)),]
  } # if(!missing(addcovar))
  if(!missing(K)) {
    K = K[match(rownames(K), colnames(expr)),
          match(colnames(K), colnames(expr))]
  } # if(!is.null(K))
  # Convert the expression data to a matrix.
  expr = as.matrix(expr)
  # Convert the SNPs into list of matrices, one for each founder.
  p2 = vector("list", dim(probs)[2])
  for(i in 1:length(p2)) {
    p2[[i]] = t(probs[,i,])
  } # for(i)
  probs = p2
  rm(p2)
  # Create an error covariance matrix.
  correctionMatrix = numeric()
  if(!missing(K)) {
    eig = eigen(K, symmetric = TRUE)
    d = eig$values
    v = eig$vectors
    if(any(d <= 0)) {
        stop("The covariance matrix is not positive definite")
    } # if(any(d <= 0))
    correctionMatrix = v %*% diag(1./sqrt(d)) %*% t(v)
    rm(eig, v, d, K)
  } # if(!missing(K))
  # Add an intercept and rotate the covariates.
  cvrt = matrix(1, nrow = 1, ncol = ncol(expr))
  if(!missing(addcovar)) {
    cvrt = rbind(cvrt, t(addcovar))
  } # if(!missing(addcovar))
  if(length(correctionMatrix) > 0) {
     cvrt = cvrt %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  q = qr(t(cvrt))
  cvrt = t(qr.Q(q))
  # Rotate and center the genes.
  if(length(correctionMatrix) > 0) {
    expr = expr %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  expr = expr - tcrossprod(expr, cvrt) %*% cvrt
  div = sqrt(rowSums(expr^2))
  div[div == 0] = 1
  expr = expr / div
  rm(div)
  # Rotate and center the SNPs. 
  for(i in 1:length(probs)) {
    if(length(correctionMatrix) > 0) {
      probs[[i]] = probs[[i]] %*% correctionMatrix
    } # if(length(correctionMatrix) > 0)
    probs[[i]] = probs[[i]] - tcrossprod(probs[[i]], cvrt) %*% cvrt
  } # for(i)
  for(i in 1:length(probs)) {
    if(i > 1) {
      for(j in 1:(i-1)) {
        probs[[i]] = probs[[i]] - rowSums(probs[[i]]*probs[[j]]) * probs[[j]]
      } # for(j)
    } # if(i > 1)
    
    div = sqrt(rowSums(probs[[i]]^2))
    drop = div < 5 * sqrt(ncol(probs[[i]]) * .Machine$double.eps)
    div[drop] = 1
    probs[[i]] = probs[[i]] / div
    probs[[i]][drop,] = 0
  } # for(i)
  probs = probs[1:(length(probs)-1)]
  rm(div, drop)
  R2 = 0;
  for(i in 1:length(probs)) {
    R2 = R2 + tcrossprod(probs[[i]], expr)^2
  } # for(i)
  print(paste("Time:", proc.time()[3] - time1, "sec."))
  # Return the LOD.
  return((-ncol(expr) * log(1.0 - R2)) / (2 * log(10)))
} # scanone.eqtl()
################################################################################
# Helper function for use in merge.analysis(). This uses the matrix eqtl 
# algorithm on a matrix of SNPs that are coded as 0, 0.5 and 1, with the 
# possibility of values in between those.
################################################################################
matrixeqtl.snps = function(pheno, geno, K, addcovar) {
  pheno = as.matrix(t(pheno))
  geno  = as.matrix(t(geno))
  # Create an error covariance matrix.
  if(!missing(K)) {
    eig = eigen(K, symmetric = TRUE)
    if(any(eig$values <= 0)) {
      stop("The covariance matrix is not positive definite")
    } # if(any(eig$values <= 0))
    correctionMatrix = eig$vectors %*% diag(1.0 / sqrt(eig$values)) %*% 
                       t(eig$vectors)
    rm(eig)
  } else {
    correctionMatrix = numeric()
  } # else
  # Add an intercept and rotate the covariates.
  cvrt = matrix(1, nrow = 1, ncol = ncol(pheno))
  if(!missing(addcovar)) {
    cvrt = rbind(cvrt, t(addcovar))
  } # if(!missing(addcovar))
  if(length(correctionMatrix) > 0) {
    cvrt = cvrt %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  q = qr(t(cvrt))
  cvrt = t(qr.Q(q))
  # Rotate and center the genes.
  if(length(correctionMatrix) > 0) {
    pheno = pheno %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  pheno = pheno - tcrossprod(pheno, cvrt) %*% cvrt
  div = sqrt(rowSums(pheno^2))
  div[div == 0] = 1
  pheno = pheno / div
  rm(div)
  # Rotate and center the SNPs. 
  if(length(correctionMatrix) > 0) {
    geno = geno %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  geno = geno - tcrossprod(geno, cvrt) %*% cvrt
  div = sqrt(rowSums(geno^2))
  drop = div < 5 * sqrt(ncol(geno) * .Machine$double.eps)
  div[drop] = 1
  geno = geno / div
  geno[drop,] = 0
  R2 = tcrossprod(geno, pheno)^2
  # Note: we return the LRS.
  return((-ncol(pheno) * log(1.0 - R2)))
} # matrixeqtl.snps()
################################################################################
# Arguments: expr: numeric matrix with samples in rows and genes in columns.
#                  Rows and columns must contain sample and gene names.
#            probs: numeric array with samples in dim[1], 8 founders in dim[2],
#                   SNPs in dim[3]. Rows, columns and slices must have names.
#            K: list containing numeric matrices with kinship coefficients. 
#                 Samples in rows and columns. Rows and columns must have names.
#            addcovar: numeric matrix with additive covariates. samples in 
#                   rows, covariates in columns. Do not include the intercept.
#                   Rows and columns must have names.
# Returns: numeric matrix with LOD values for each gene and snp.
scanone.loo = function(expr, probs, K, addcovar, snps, sex) {
  # expr, probs and snps cannot be null.
  if(missing(expr) || missing(probs) || missing(K)) {
    stop(paste("One of the required arguments is null.  Please make sure",
         "that expr, probs and K are not null."))
  } # if(is.null(expr) | ...
  expr = t(expr)
  # Verify that the number of rows and columns match in the data.
  # Same number of samples in expr, probs and addcovar.
  if(ncol(expr) != dim(probs)[[1]]) {
    stop(paste("The number of columns in the expression matrix (",
         ncol(expr),") is not equal to the number of rows in probs (",
         dim(probs)[[1]],")."))
  } # if(ncol(expr) != dim(probs)[[1]])
  
  # Verify that the number of samples is the same in addcovar and expr.
  if(!missing(addcovar)) {
    if(nrow(addcovar) != ncol(expr)) {
      stop(paste("The number of rows in addcovar (", nrow(addcovar),") must",
           "equal the number of columns in expr (", nrow(expr) ,")."))
    } # if(nrow(addcovar) != nrow(expr))
  } # if(!missing(addcovar))
  # Verify that the number of samples is the same in the kinship and 
  # expression matrices.
  if(!is.list(K)) stop("K must be a list of matrices in scanone.loo().")
  if(!all(sapply(K, nrow) == ncol(expr))) {
    paste(stop("The number of rows in the kinship matrix (", nrow(K)[[1]],
          ") does not equal the number of columns in expr (", ncol(expr),
          ")."))
  } # if(nrow(K) != ncol(expr))
  if(!all(sapply(K, nrow) == ncol(expr))) {
    paste(stop("The number of columns in the kinship matrix (", ncol(K)[[1]],
          ") does not equal the number of columns in expr (", ncol(expr),
           ")."))
  } # if(ncol(K) != ncol(expr))
  # Match up the sample IDs.
  if(is.null(colnames(expr))) {
    stop(paste("expr must have sample IDs in the column names."))
  } # if(is.null(colnames(expr)))
  probs = probs[match(dimnames(probs)[[1]], colnames(expr)),,]
  if(!missing(addcovar)) {
    addcovar = addcovar[match(rownames(addcovar), colnames(expr)),]
  } # if(!missing(addcovar))
  for(i in 1:length(K)) {
    K[[i]] = K[[i]][match(rownames(K[[i]]), colnames(expr)),
             match(colnames(K[[i]]), colnames(expr))]
  } # for(i)
  # Convert the expression data to a matrix.
  expr = as.matrix(expr)
  # Convert the SNPs into list of matrices, one for each founder.
  p2 = vector("list", dim(probs)[2])
  for(i in 1:length(p2)) {
    p2[[i]] = t(probs[,i,])
  } # for(i)
  probs = p2
  rm(p2)
  R2 = rep(0, nrow(snps))
  for(c in 1:length(K)) {
    # Get the SNP range for the current chormosome.
    snprng = which(snps[,2] == names(K)[c])
    # Create error covariance matrices.
    correctionMatrix = numeric()
    eig = eigen(K[[c]], symmetric = TRUE)
    d = eig$values
    v = eig$vectors
    if(any(d <= 0)) {
       stop("The covariance matrix is not positive definite")
    } # if(any(d <= 0))
    correctionMatrix = v %*% diag(1.0 / sqrt(d)) %*% t(v)
    rm(eig, d, v)
    # Add an intercept and rotate the covariates.
    cvrt = matrix(1, nrow = 1, ncol = ncol(expr))
    if(!missing(addcovar)) {
      cvrt = rbind(cvrt, t(addcovar))
    } # if(!missing(addcovar))
    cvrt = cvrt %*% correctionMatrix
    q    = qr(t(cvrt))
    cvrt = t(qr.Q(q))
    # Rotate and center the genes. e will be the chromosome specific expr.
    e = expr %*% correctionMatrix
    e = e - tcrossprod(e, cvrt) %*% cvrt
    div = sqrt(rowSums(e^2))
    div[div == 0] = 1
    e = e / div
    rm(div)
    # Rotate and center the SNPs. p will be the chromosome specific probs.
    p = vector("list", length(probs))
    for(i in 1:length(probs)) {
      p[[i]] = probs[[i]][snprng,] %*% correctionMatrix
      p[[i]] = p[[i]] - tcrossprod(p[[i]], cvrt) %*% cvrt
    } # for(i)
    for(i in 1:length(p)) {
      if(i > 1) {
        for(j in 1:(i-1)) {
          p[[i]] = p[[i]] - rowSums(p[[i]] * p[[j]]) * p[[j]]
        } # for(j)
      } # if(i > 1)
    
      div = sqrt(rowSums(p[[i]]^2))
      drop = div < 5 * sqrt(ncol(p[[i]]) * .Machine$double.eps)
      div[drop] = 1
      p[[i]] = p[[i]] / div
      p[[i]][drop,] = 0
    } # for(i)
    p = p[1:(length(p)-1)]
    rm(div, drop)
    for(i in 1:length(p)) {
      R2[snprng] = R2[snprng] + tcrossprod(p[[i]], e)^2
    } # for(i)
  } # for(c)
  # Return the LOD.
  return((-ncol(expr) * log(1.0 - R2)) / (2 * log(10)))
} # scanone.loo()
