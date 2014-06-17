################################################################################
# Association mapping using Sanger SNPs imputed onto DO genomes using the 
# DO haplotype contributions. Kinship is scaled and used as in QTLRel.
# Daniel Gatti
# Dan.Gatti@jax.org
# Jan. 18, 2012
################################################################################
# Arguments: pheno: phenotype, data.frame with phenotypes.
#            pheno.col: numeric or character. If numeric, the column number in 
#                       the phenotype matrix to use. If character, the column
#                       name in the phenotype matrix to use.
#            probs: 3D numeric array, with founder probabilities for all 
#                   samples.
#            K: numeric matrix, with kinship for samples.
#            addcovar: numeric matrix, with fixed covariates.
#            snps: data.frame, with MUGA SNPs.
#            chr: character, the chromosome to use.
#            start: numeric, the start location along the chromosome. Values
#                   greater than 200 will be interpreted in bp. Values less
#                   than or equal to 200 will be interpreted as MB.
#            end: numeric, the end location along the chromosome. Values
#                 greater than 200 will be interpreted in bp. Values less
#                 than or equal to 200 will be interpreted as MB.
#            model: character, one of "additive", "dominance" or "full" that is 
#                   the model to fit.
#            scan: character, one of "one" or "two" that indicates whether to
#                  fit a single locus or all pairs of loci.
#            snp.file: character, full path to a file of SNP, Indels, structural
#                      variants (or a combination of these) zipped and tabix
#                      indexed.
assoc.map = function(pheno, pheno.col = 1, probs, K, addcovar, snps,
                 chr, start, end, model = c("additive", "dominance", "full"),
                 scan = c("one", "two"), output = c("lod", "p-value", "bic"),
                 snp.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz") {
  scan  = match.arg(scan)
  model = match.arg(model)
  output = match.arg(output)
  if(scan == "two") {
    stop("merge.analysis: scan two not implemented yet.")
  } # if(scan == "two")
  if(model != "additive") {
    stop("merge.analysis: Dominance and full models not implemented yet.")
  } # if(model != "additive")
  if(missing(pheno)) {
    stop("The phenotypes are missing. Please supply a phenotype data frame.")
  } # if(missing(pheno))
  if(missing(probs)) {
    stop(paste("The founder probabilities are missing. Please supply an array",
         "of founder probabilities."))
  } # if(missing(probs))
  if(missing(K)) {
    stop(paste("The kinship matrix is missing. Please supply a kinship matrix",
         "for all samples."))
  } # if(missing(K))
  if(missing(chr)) {
    stop("The QTL chromosome must be specified.")
  } # if(missing(chr))
  if(missing(start)) {
    stop("The QTL start position must be specified.")
  } # if(missing(start))
  if(missing(end)) {
    stop("The QTL end position must be specified.")
  } # if(missing(end))
  # Subset all of the data to only contain the samples in common.
  pheno = pheno[!is.na(pheno[,pheno.col]),,drop = FALSE]
  samples = intersect(rownames(pheno), dimnames(probs)[[1]])
  if(length(samples) == 0) {
    stop(paste("rownames(pheno) does not contain any sample IDs that",
         "match dimnames(probs)[[1]]. Please verify that both variables",
         "contain the same sample IDs."))
  } # if(length(samples) == 0)
  if(!missing(addcovar)) {
    addcovar = addcovar[rownames(addcovar) %in% samples,,drop = FALSE]
    if(nrow(addcovar) == 0) {
      stop(paste("There are no common sample IDs between rownames(addcovar)",
           "and rownames(pheno). Please ensure that there is some overlap",
           "in sample IDs between rownames(addcovar) and rownames(pheno)."))
    } # if(nrow(addcovar) == 0)
    remove = which(rowSums(is.na(addcovar)) > 0)
    if(length(remove) > 0) {
      samples = samples[-remove]
      addcovar = addcovar[-remove,,drop = FALSE]
    } # if(length(remove) > 0)
  } # if(!missing(addcovar))
  pheno = pheno[rownames(pheno) %in% samples,,drop = FALSE]
  probs = probs[dimnames(probs)[[1]] %in% samples,,]
  probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,]
  K = K[rownames(K) %in% samples, colnames(K) %in% samples]
  K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), colnames(K))]
  print(paste("Mapping with", nrow(pheno), "samples."))
  stopifnot(all(rownames(pheno) == dimnames(probs)[[1]]))
  stopifnot(all(rownames(pheno) == rownames(K)))
  stopifnot(all(rownames(pheno) == colnames(K)))
  probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
  if(!all(snps[,1] == dimnames(probs)[[3]])) {
    stop(paste("All of the SNP IDs in snps do not match the SNP IDs in",
         "dimnames(probs)[[3]]."))
  } # if(!all(snps[,1] == dimnames(probs)[[3]]))
  # Subset the SNPs to be the same.
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
  stopifnot(all(snps[,1] == dimnames(probs)[[3]]))
  start = as.numeric(start)
  if(start <= 200) {
    start = start * 1e6
  } # if(start <= 200)
  end = as.numeric(end)
  if(end <= 200) {
    end = end * 1e6
  } # if(end <= 200)
  # Expand the start and end to the nearest markers.
  snps = snps[snps[,2] == chr,]
  snps = snps[intersect(which(snps[,3] >= start * 1e-6) - 1,
                        which(snps[,3] <= end   * 1e-6) + 1),]
  probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  
  gr = GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  print("Retrieving SNPs...")
  con = TabixFile(snp.file)
  open(con)
  # Get the column names from the last row of the header info.
  hdr = headerTabix(con)
  hdr = strsplit(hdr$header, split = "\t")[[length(hdr$header)]]
  hdr = sub("^#", "", hdr)
  # Retrieve the data.
  sanger = scanTabix(con, param = gr)
  # Close the connection.
  close(con)
  print(paste("Retrieved", length(sanger[[1]]), "SNPs."))
  # Split the columns up (tab-delimited).
  print("Finding unique SNP patterns...")
  sanger = strsplit(sanger[[1]], split = "\t")
  sanger = matrix(unlist(sanger), length(sanger), length(sanger[[1]]),
           dimnames = list(sanger$POS, hdr), byrow = TRUE)
  
  # Keep only high quality reads where all quality scores = 1.
  qual.columns = grep("quality", colnames(sanger))
  sanger = sanger[rowMeans(sanger[,qual.columns] == "1") == 1,,drop = FALSE]
  # Remove quality scores.
  sanger = sanger[,-qual.columns]
  # Remove SNPs with "NN" calls because we don't know how to assign
  # the founder alleles in that case.
  sanger = sanger[rowSums(sanger == "NN") == 0,]
  # Make reference allele.
  ref.col = which(colnames(sanger) == "REF")
  ref = paste(sanger[,ref.col], sanger[,ref.col], sep = "")
  pos = sanger[,colnames(sanger) == "POS"]
  # Convert allele calls to numbers.
  sanger = matrix(as.numeric(sanger[,-1:-5] != ref), nrow(sanger), ncol(sanger) - 5,
           dimnames = list(NULL, colnames(sanger)[-1:-5]))
  # Remove SNPs that are all 0s.
  na.rows = which(rowSums(is.na(sanger)) > 0)
  remove = na.rows[rowSums(sanger[na.rows,], na.rm = TRUE) == 0]
  if(length(remove) > 0) {
    sanger = sanger[-remove,]
    ref = ref[-remove]
    pos = pos[-remove]
  } # if(length(remove) > 0)
  # Change the allele calls so that the MAF <= 4.
  wh = which(rowSums(sanger) > 4)
  sanger[wh,] = 1 - sanger[wh,]
  # Reorder the strains to match DOQTL founders.
  m = match(c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ",
              "CAST/EiJ", "PWK/PhJ", "WSB/EiJ"), colnames(sanger))
  if(any(is.na(m))) {
    stop("One of the DO founders is not in the SNP file.")
  } # if(any(is.na(m)))
  sanger = sanger[,m]
  rownames(sanger) = pos
  # Get the SDPs (this is DO/CC specific).
  # Get the SDPs by treating the SNP patterns as binary numbers.
  vec = 2^((ncol(sanger)-1):0)
  sdps = tcrossprod(vec, sanger)[1,]
  sdps = data.frame(pos = as.numeric(pos), sdps, sanger)
  print("Calculating mapping statistic")
  retval = NULL
  if(scan == "one") {
    retval = assoc.scan1(pheno = pheno, pheno.col = pheno.col, probs = probs, 
             K = K, addcovar = addcovar, sdps = sdps, snps = snps, 
             model = model, output = output)
  } else {
# TBD: implement pair scan.
    retval = assoc.scan2(pheno, pheno.col, addcovar, sanger, pos, sdps, snps)
  } # else
  colnames(retval) = sub("^sdp\\.", "", colnames(retval))
  attr(retval, "class") = c(attr(retval, "class"), "assoc")
  return(retval)
} # assoc.map()
##########
# Association mapping permutations.
assoc.map.perms = function(pheno, pheno.col = 1, probs, addcovar, snps,
            model = c("additive", "dominance", "full"),
            scan = c("one", "two"), output = c("lod", "p-value", "bic"),
            snp.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz",
            nperm = 1000) {
  scan  = match.arg(scan)
  model = match.arg(model)
  output = match.arg(output)
  if(scan == "two") {
    stop("merge.analysis: scan two not implemented yet.")
  } # if(scan == "two")
  if(model != "additive") {
    stop("merge.analysis: Dominance and full models not implemented yet.")
  } # if(model != "additive")
  if(missing(pheno)) {
    stop("The phenotypes are missing. Please supply a phenotype data frame.")
  } # if(missing(pheno))
  if(missing(probs)) {
    stop(paste("The founder probabilities are missing. Please supply an array",
         "of founder probabilities."))
  } # if(missing(probs))
  # Subset all of the data to only contain the samples in common.
  pheno = pheno[!is.na(pheno[,pheno.col]),,drop = FALSE]
  samples = intersect(rownames(pheno), dimnames(probs)[[1]])
  if(length(samples) == 0) {
    stop(paste("rownames(pheno) does not contain any sample IDs that",
         "match dimnames(probs)[[1]]. Please verify that both variables",
         "contain the same sample IDs."))
  } # if(length(samples) == 0)
  if(!missing(addcovar)) {
    addcovar = addcovar[rownames(addcovar) %in% samples,,drop = FALSE]
    if(nrow(addcovar) == 0) {
      stop(paste("There are no common sample IDs between rownames(addcovar)",
           "and rownames(pheno). Please ensure that there is some overlap",
           "in sample IDs between rownames(addcovar) and rownames(pheno)."))
    } # if(nrow(addcovar) == 0)
    remove = which(rowSums(is.na(addcovar)) > 0)
    if(length(remove) > 0) {
      samples = samples[-remove]
      addcovar = addcovar[-remove,,drop = FALSE]
    } # if(length(remove) > 0)
  } # if(!missing(addcovar))
  pheno = pheno[rownames(pheno) %in% samples,,drop = FALSE]
  probs = probs[dimnames(probs)[[1]] %in% samples,,]
  probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,]
  print(paste("Mapping with", nrow(pheno), "samples."))
  stopifnot(all(rownames(pheno) == dimnames(probs)[[1]]))
  # Subset the SNPs to be the same.
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
  stopifnot(all(snps[,1] == dimnames(probs)[[3]]))
  chr = unique(snps[,2])
  start = 0
  end = 200e6
  retval = matrix(0, length(chr), nperm)
  # Loop through each chromosome.
  for(c in 1:length(chr)) {
    # Expand the start and end to the nearest markers.
    curr.snps = snps[snps[,2] == chr[c],]
    snps.to.keep = which(dimnames(probs)[[3]] %in% curr.snps[,1])
  
    gr = GRanges(seqnames = chr[c], ranges = IRanges(start = start, end = end))
    print(paste("Retrieving SNPs on Chr", chr[c], "..."))
    con = TabixFile(snp.file)
    open(con)
    # Get the column names from the last row of the header info.
    hdr = headerTabix(con)
    hdr = strsplit(hdr$header, split = "\t")[[length(hdr$header)]]
    hdr = sub("^#", "", hdr)
    # Retrieve the data.
    sanger = scanTabix(con, param = gr)
    # Close the connection.
    close(con)
    print(paste("Retrieved", length(sanger[[1]]), "SNPs."))
    # Split the columns up (tab-delimited).
    print("Finding unique SNP patterns...")
    sanger = strsplit(sanger[[1]], split = "\t")
    sanger = matrix(unlist(sanger), length(sanger), length(sanger[[1]]),
             dimnames = list(sanger$POS, hdr), byrow = TRUE)
  
    # Keep only high quality reads where all quality scores = 1.
    qual.columns = grep("quality", colnames(sanger))
    sanger = sanger[rowMeans(sanger[,qual.columns] == "1") == 1,,drop = FALSE]
    # Remove quality scores.
    sanger = sanger[,-qual.columns]
    # Remove SNPs with "NN" calls because we don't know how to assign
    # the founder alleles in that case.
    sanger = sanger[rowSums(sanger == "NN") == 0,]
    # Make reference allele.
    ref.col = which(colnames(sanger) == "REF")
    ref = paste(sanger[,ref.col], sanger[,ref.col], sep = "")
    pos = sanger[,colnames(sanger) == "POS"]
    # Convert allele calls to numbers.
    sanger = matrix(as.numeric(sanger[,-1:-5] != ref), nrow(sanger), ncol(sanger) - 5,
             dimnames = list(NULL, colnames(sanger)[-1:-5]))
    # Remove SNPs that are all 0s.
    na.rows = which(rowSums(is.na(sanger)) > 0)
    remove = na.rows[rowSums(sanger[na.rows,], na.rm = TRUE) == 0]
    if(length(remove) > 0) {
      sanger = sanger[-remove,]
      ref = ref[-remove]
      pos = pos[-remove]
    } # if(length(remove) > 0)
    # Change the allele calls so that the MAF <= 4.
    wh = which(rowSums(sanger) > 4)
    sanger[wh,] = 1 - sanger[wh,]
    # Reorder the strains to match DOQTL founders.
    m = match(c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ",
                "CAST/EiJ", "PWK/PhJ", "WSB/EiJ"), colnames(sanger))
    if(any(is.na(m))) {
      stop("One of the dO founders is not in the SNP file.")
    } # if(any(is.na(m)))
    sanger = sanger[,m]
    rownames(sanger) = pos
    # Get the SDPs (this is DO/CC specific).
    # Get the SDPs by treating the SNP patterns as binary numbers.
    vec = 2^((ncol(sanger)-1):0)
    sdps = tcrossprod(vec, sanger)[1,]
    sdps = data.frame(pos = as.numeric(pos), sdps, sanger)
    print("Calculating mapping statistic")
    if(scan == "one") {
      retval[c,] = assoc.scan1.perms(pheno = pheno, pheno.col = pheno.col,
               probs = probs[,,snps.to.keep], addcovar = addcovar,
               sdps = sdps, snps = curr.snps, model = model, output = output,
               nperm = nperm)
    } else {
# TBD: implement pair scan.
      retval[c,] = assoc.scan2(pheno, pheno.col, addcovar, sanger, pos, sdps, snps)
    } # else
    
  } # for(c)
  if(output == "p-value") {
    retval = apply(retval, 2, min)
  } else {
    retval = apply(retval, 2, max)
  } # else
  
  attr(retval, "class") = c(attr(retval, "class"), "assoc")
  return(retval)
} # assoc.map.perms()
###
# Helper function for one-way scan.
# The chromosome range is based on the start and end of the Sanger SNP locations
# passed in by the user.
assoc.scan1 = function(pheno, pheno.col, probs, K, addcovar, sdps, 
              snps, model, output) {
  # Get the error covariance matrix from QTLRel.
  err.cov = diag(nrow(pheno))
  if(!missing(K)) {
    vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
                EE = diag(nrow(pheno)))
    vc = estVC(y = pheno[,pheno.col], x = addcovar, v = vTmp)
    err.cov = matrix(0, nrow(K), ncol(K))
    for(j in which(vc$nnl)) {
      err.cov = err.cov + vTmp[[j]] * vc$par[names(vTmp)[j]]
    } # for(j)
    rm(vTmp)
  } # if(!missing(K))
  BIC.full = NULL
  if(output == "bic") {
    # Get the SS.full (8 founder state).
    tmp = probs[,,-dim(probs)[3]]
    for(i in 1:dim(tmp)[3]) {
      tmp[,,i] = 0.5 * (probs[,,i] + probs[,,i+1])
    } # for(i)
    qtl.full = scanone(pheno = pheno, pheno.col = pheno.col, 
               probs = tmp, K = K, addcovar = addcovar, 
               snps = snps[-nrow(snps),])
    rm(tmp)
    BIC.full = -2 * qtl.full$lod$A$lod + (ncol(addcovar) +  dim(probs)[[2]]) * 
               log(nrow(pheno))
  } # if(output == "bic")
  # Convert the DO haplotype probabilities to Sanger alleles.
  unique.geno = dohap2sanger(probs, snps, sdps)
  keep = which(!is.na(pheno[,pheno.col]))
  if(!missing(addcovar)) {
    keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
    if(!missing(K)) {
      lrs = matrixeqtl.snps(pheno = pheno[keep,pheno.col], 
            geno = unique.geno$allele.probs[keep,], 
            K = err.cov[keep,keep], addcovar = addcovar[keep,])
    } else {
      lrs = matrixeqtl.snps(pheno = pheno[keep,pheno.col],
            geno = unique.geno$allele.probs[keep,], addcovar = addcovar[keep,])
    } # else
  } else {
    if(!missing(K)) {
      lrs = matrixeqtl.snps(pheno = pheno[keep,pheno.col],
            geno = unique.geno$allele.probs[keep,], K = err.cov[keep,keep])
    } else {
      lrs = matrixeqtl.snps(pheno = pheno[keep,pheno.col],
            geno = unique.geno$allele.probs[keep,])
    } # else
  } # else
  if(output == "bic") {
    BIC.reduced = -2 * lrs + (ncol(addcovar) + 3) * log(nrow(pheno))
    stat = BIC.reduced[unique.geno$map]
    bfull = rep(0, nrow(sdps))
    for(i in 1:(nrow(snps) - 1)) {
      bfull[sdps[,1] >= snps[i,3] * 1e6 & sdps[,1] <= snps[i+1,3] * 1e6] = BIC.full[i]
    } # for(i)
    return(data.frame(chr = rep(snps[1,2], length(stat)), 
           sdp = sdps, BIC.full = bfull, BIC.red = stat))
  } else if(output == "p-value") {
    stat = pchisq(q = lrs, df = 1, lower.tail = FALSE)[unique.geno$map]
    return(data.frame(chr = rep(snps[1,2], length(stat)), 
           sdp = sdps, p.value = stat))
  } else if(output == "lod") {
    stat = (lrs / (2 * log(10)))[unique.geno$map]
    return(data.frame(chr = rep(snps[1,2], length(stat)), 
           sdp = sdps, LOD = stat))
  } # else if(output == "lod")
  
} # assoc.scan1
##########
# Scan one permutations. 
assoc.scan1.perms = function(pheno, pheno.col, probs, K, addcovar,
                    sdps, snps, model, output, nperm = 1000) {
  BIC.full = NULL
  if(output == "bic") {
    # Get the SS.full (8 founder state).
    tmp = probs[,,-dim(probs)[3]]
    for(i in 1:dim(tmp)[3]) {
      tmp[,,i] = 0.5 * (probs[,,i] + probs[,,i+1])
    } # for(i)
    qtl.full = scanone(pheno = pheno, pheno.col = pheno.col, 
               probs = tmp, K = K, addcovar = addcovar, 
               snps = snps[-nrow(snps),])
    rm(tmp)
    BIC.full = -2 * qtl.full$lod$A$lod + (ncol(addcovar) +  dim(probs)[[2]]) * 
               log(nrow(pheno))
  } # if(output == "bic")
  # Convert the DO haplotype probabilities to Sanger alleles.
  unique.geno = dohap2sanger(probs, snps, sdps)
  keep = which(!is.na(pheno[,pheno.col]))
  if(!missing(addcovar)) {
    # If we have additive covariates, we must permute them with the phenotype.
    keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
    ph = pheno[keep,pheno.col]
    cv = addcovar[keep,,drop=FALSE]
    max.lrs = rep(0, nperm)
    perms = matrix(1:length(ph), nrow = length(ph), ncol = nperm)
    perms = apply(perms, 2, sample)
    for(i in 1:nperm) {
      lrs = matrixeqtl.snps(pheno = ph[perms[,i]],
            geno = unique.geno$allele.probs[keep,], addcovar = cv[perms[,i],,drop=FALSE])
      max.lrs[i] = max(lrs[,1])
    } # for(i)
  } else {
    # In this case, we can run all of the permutations at once since there is
    # no additive covariate to permute.
    ph = pheno[keep,pheno.col]
    perm.pheno = matrix(ph, nrow = length(ph), ncol = nperm)
    perm.pheno = apply(perm.pheno, 2, sample)
    lrs = matrixeqtl.snps(pheno = perm.pheno,
          geno = unique.geno$allele.probs[keep,,drop=FALSE])
    max.lrs = apply(lrs, 2, max)
  } # else
  if(output == "bic") {
stop("Permutations for BIC not implemented yet.")
# TBD: We have to permuate the full haplotype model as well.
    BIC.reduced = -2 * lrs + (ncol(addcovar) + 3) * log(nrow(pheno))
    stat = BIC.reduced[unique.geno$map]
    bfull = rep(0, nrow(sdps))
    for(i in 1:(nrow(snps) - 1)) {
      bfull[sdps[,1] >= snps[i,3] * 1e6 & sdps[,1] <= snps[i+1,3] * 1e6] = BIC.full[i]
    } # for(i)
    return(data.frame(chr = rep(snps[1,2], length(stat)), 
           sdp = sdps, BIC.full = bfull, BIC.red = stat))
  } else if(output == "p-value") {
    stat = pchisq(q = max.lrs, df = 1, lower.tail = FALSE)
    return(stat)
  } else if(output == "lod") {
    stat = max.lrs / (2 * log(10))
    return(stat)
  } # else if(output == "lod")
    
} # assoc.scan1.perms()
###
# Helper function for two-way scan.
assoc.scan2 = function(pheno, pheno.col, probs, K, addcovar, sdps, 
              snps, model, output) {
  stop("merge.scan2 not implemented yet.")
  # Get the error covariance matrix from QTLRel.
  err.cov = diag(nrow(pheno))
  if(!missing(K)) {
    vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
                EE = diag(nrow(pheno)))
    vc = NULL
    if(missing(addcovar)) {
      vc = estVC(y = pheno[,pheno.col], v = vTmp)
    } else {
      vc = estVC(y = pheno[,pheno.col], x = addcovar, v = vTmp)
    } # else
    err.cov = matrix(0, nrow(K), ncol(K))
    for(j in which(vc$nnl)) {
      err.cov = err.cov + vTmp[[j]] * vc$par[names(vTmp)[j]]
    } # for(j)
    rm(vTmp)
  } # if(!missing(K))
  # Convert the DO haplotype probabilities to Sanger alleles.
  unique.geno = dohap2sanger(probs, snps, sdps)
  keep = which(!is.na(pheno[,pheno.col]))
  if(!missing(addcovar)) {
    keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
    if(!missing(K)) {
      lrs = matrix(0, ncol(unique.geno$allele.probs), ncol(unique.geno$allele.probs))
      for(i in 1:ncol(unique.geno$allele.probs)) {
        addcovar = cbind(addcovar, unique.geno$allele.probs[,i])
        lrs[i,] = matrixeqtl.snps(pheno = pheno[keep,pheno.col], 
                  geno = unique.geno$allele.probs[keep,], 
                  K = err.cov[keep,keep], addcovar = addcovar[keep,,drop=FALSE])
      } # for(i)
    } else {
      lrs = matrix(0, ncol(unique.geno$allele.probs), ncol(unique.geno$allele.probs))
      for(i in 1:ncol(unique.geno$allele.probs)) {
        addcovar = cbind(addcovar, unique.geno$allele.probs[,i])
        lrs[i,] = matrixeqtl.snps(pheno = pheno[keep,pheno.col],
                  geno = unique.geno$allele.probs[keep,], addcovar = 
                  addcovar[keep,,drop=FALSE])
      } # for(i)
    } # else
  } else {
    if(!missing(K)) {
      lrs = matrix(0, ncol(unique.geno$allele.probs), ncol(unique.geno$allele.probs))
      for(i in 1:ncol(unique.geno$allele.probs)) {
        addcovar = as.matrix(unique.geno$allele.probs[,i])
        lrs[i,] = matrixeqtl.snps(pheno = pheno[keep,pheno.col],
                  geno = unique.geno$allele.probs[keep,], K = err.cov[keep,keep],
                  addcovar = addcovar[keep,,drop=FALSE])
      } # for(i)
    } else {
      lrs = matrix(0, ncol(unique.geno$allele.probs), ncol(unique.geno$allele.probs))
      for(i in 1:ncol(unique.geno$allele.probs)) {
        addcovar = as.matrix(unique.geno$allele.probs[,i])
        lrs[i,] = matrixeqtl.snps(pheno = pheno[keep,pheno.col],
                  geno = unique.geno$allele.probs[keep,], addcovar = 
                  addcovar[keep,,drop=FALSE])
      } # for(i)
    } # else
  } # else
  if(output == "lod") {
stop("Return a matrix of LOD scores")
    stat = (lrs / (2 * log(10)))[unique.geno$map]
    return(data.frame(chr = rep(snps[1,2], length(stat)), 
           sdp = sdps, LOD = stat))
  } 
} # assoc.scan2
###
# Map Sanger SNPs onto DO genomes, given the 8 founder haplotyep contributions.
# Arguments: probs: 3D array of founder haplotype contributions. num.samples x
#                   num.founders x num.markers.
#            snps: data.frame containing genotyping array marker locations. 
#                  SNP ID, chr, Mb and cM in columns 1:4.
#            sdps: data.frame containing bp, SDP and Sanger SNPs.
dohap2sanger = function(probs, snps, sdps) {
  # Get mean haplotype contributions between each pair of markers.
  mean.haps = sapply(1:(dim(probs)[3] - 1), function(z) { 
                0.5 * (probs[,,z] + probs[,,z+1])
              })
  mean.haps = array(mean.haps, c(dim(probs)[1], dim(probs)[2], dim(probs)[3] - 1))
  # Create bp position, SDPs and a numeric matrix of sanger SNPs.
  pos = sdps[,1]
  sanger = as.matrix(sdps[,-1:-2])
  sdps = sdps[,2]
  snps[,3] = snps[,3] * 1e6
  # Function to take the Sanger SNPs, a range, and return unique SDPs
  # and a map from the Sanger SNPs to the SDPs.
  map.fxn = function(start, end, haps) {
              rng = which(pos >= start & pos < end)
              map = NULL
              allele.probs = NULL
              if(length(rng) > 0) {
                sdp.rng = which(!duplicated(sdps[rng]))
                unique.sdp = sdps[rng][sdp.rng]
                unique.sanger = sanger[rng,,drop = FALSE][sdp.rng,,drop = FALSE]
                map = match(sdps[rng], unique.sdp)
                allele.probs = tcrossprod(haps, unique.sanger)                
              } # if(length(rng) > 0)
              return(list(map = map, allele.probs = allele.probs))
            } # function()
  # Create space for the unique SDPs and alleles that is 10 times smaller
  # than the total number of Sanger SNPs.
  map = rep(NA, length(sdps))
  allele.probs = matrix(0, dim(probs)[1], length(sdps) * 0.1, dimnames =
                 list(dimnames(probs)[[1]], NULL))
  map.index   = 1
  probs.index = 1
  # Process positions before the first marker.
  if(pos[1] < snps[1,3]) {
    tmp = map.fxn(pos[1], snps[1,3], probs[,,1])
    if(length(tmp$map) > 0) {
      map[map.index:(map.index + length(tmp$map) - 1)] = tmp$map + probs.index - 1
      allele.probs[,probs.index:(probs.index + ncol(tmp$allele.probs) - 1)] = tmp$allele.probs
      map.index = map.index + length(tmp$map)
      probs.index = probs.index + ncol(tmp$allele.probs)
    } # if(length(tmp$map) > 0)
  } # if(pos[1] < snps[1,3])
  # Process positions between markers.
  for(i in 1:(nrow(snps)-1)) {
    tmp = map.fxn(snps[i,3], snps[i+1,3], mean.haps[,,i])
    if(length(tmp$map) > 0) {
      map[map.index:(map.index + length(tmp$map) - 1)] = tmp$map + probs.index - 1
      allele.probs[,probs.index:(probs.index + ncol(tmp$allele.probs) - 1)] = tmp$allele.probs
      map.index = map.index + length(tmp$map)
      probs.index = probs.index + ncol(tmp$allele.probs)
    } # if(length(tmp$map) > 0)
  } # for(i)
  # Process positions past the last marker.
  if(pos[length(pos)] > snps[nrow(snps),3]) {
    tmp = map.fxn(snps[nrow(snps),3], pos[length(pos)]+1, probs[,,dim(probs)[3]])
    if(length(tmp$map) > 0) {
      map[map.index:(map.index + length(tmp$map) - 1)] = tmp$map + probs.index - 1
      allele.probs[,probs.index:(probs.index + ncol(tmp$allele.probs) - 1)] = tmp$allele.probs
      map.index = map.index + length(tmp$map)
      probs.index = probs.index + ncol(tmp$allele.probs)
    } # if(length(tmp$map) > 0)
  } # if(pos[1] < snps[1,3])
  allele.probs = allele.probs[,1:probs.index]
  return(list(map = map, allele.probs = allele.probs))
} # dohap2sanger()
###
# Plot for one-way merge analysis scan.
# Arguments: results: data.frame from assoc.map.
#            mgi.file: character string containing the full path to the MGI
#                      feature file.
#            highlight: character vector of gene symbols to highlight.
#            highlight.col: color vector of highlight colors.
#            thr: LOD score to color red in the plot. All SNPs with a LOD
#                 above this threshold will be returned.
#            ...: arguments to be passed to plot.
assoc.plot = function(results, 
  mgi.file = "ftp://ftp.jax.org/SNPtools/genes/MGI.20130703.sorted.txt.gz",
  highlight, highlight.col = "red", thr, ...) {
  old.par = par(no.readonly = TRUE)
  start = results$pos[1] * 1e-6
  end   = results$pos[nrow(results)] * 1e-6
  call = match.call()
  if("xlim" %in% names(call)) {
     xlim = eval(call$xlim)
print(xlim)
     results = results[results[,2] >= xlim[1] * 1e6 & results[,2] <= xlim[2] * 1e6,]
     if(any(xlim > 200)) {
       xlim = xlim * 1e-6
     } # if(xlim > 200)
     call$xlim = xlim
  } # if("xlim" %in% names(call))
  type = "bic"
  if(colnames(results)[ncol(results)] == "LOD") {
    diff = results$LOD
    type = "lod"
  } else if(colnames(results)[ncol(results)] == "p.value") {
    diff = -log10(results$p.value)
    type = "pv"
  } else {
    diff = results$BIC.full - results$BIC.red
  } # else
  # Color the points with statistic > thr in red.
  col = rep(1, nrow(results))
  if(!missing(thr)) {
    points.ge.thr = which(diff >= thr)
    col[points.ge.thr] = 2
  } # else
  layout(matrix(1:2, 2, 1), heights = c(0.4, 0.6))
  par(plt = c(0.12, 0.99, 0, 0.9), las = 1)
  plot(results$pos * 1e-6, diff, ann = FALSE, xaxt = "n", pch = 20, col = col, ...)
  usr = par("usr")
  par(las = 3)
  txt = "BIC.full - BIC.reduced"
  if(type == "pv") {
    txt = "-log10(p-value)"
  } else if(type == "lod") {
    txt = "LOD"
  } # else
  mtext(side = 2, line = 3, text = txt)
  par(las = 1)
  par(plt = c(0.12, 0.99, 0.15, 1))
  mgi = get.mgi.features(file = mgi.file, chr = results[1,1],
        start = start, end = end, source = "MGI",
        type = "gene")
  col = "black"
  if(!missing(highlight)) {
    col = rep("black", nrow(mgi))
    m = match(highlight, mgi$Name)
    col[m] = highlight.col
  } # if(!missing(highlight))
  genes = gene.plot(mgi = mgi, col = col, xlim = usr[1:2])
  if(!missing(thr)) {
    return(results[points.ge.thr,])
  } else {
    return(results)
  } # else
  par(old.par)
} # assoc.plot()
# Plot specific SDP locations along a chromosome.
# This idea is from Yalcin et.al., Genetics, 2005.
# results: the data.frame output from merge.analysis().
# sdps: strings of 0s and 1s that match some of the SDPs in column
#       2 of the merge.analysis() results.
# ...: additional arguments to pass to plot.
sdp.plot = function(results, sdps, ...) {
  plot(results[,1], rep(1, nrow(results)), col = 0, ylim = c(0, length(sdps)),
       yaxs = "i", ann = FALSE, axes = FALSE, ...)
  abline(h = 0:length(sdps))
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  mtext(side = 2, at = 1:length(sdps) - 0.5, text = sdps, line = 0.1, 
        adj = 1, las = 2)
  ax = axis(side = 1, labels = FALSE)
  axis(side = 1, at = ax, labels = ax * 1e-6)
  # Find the SDPs in the results.
  for(i in 1:length(sdps)) {
    m = results[which(results[,2] == sdps[i]),1]
    y = matrix(rep((i-1):(i), length(m)), nrow = 2)
    y = rbind(y, rep(NA, ncol(y)))
    y = as.vector(y)
    x = matrix(rep(m, each = 2), nrow = 2)
    x = rbind(x, rep(NA, ncol(x)))
    x = as.vector(x)
    lines(x, y)
  } # for(i)
} # sdp.plot()
