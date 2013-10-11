################################################################################
# Merge analysis as described in Yalcin et.al., Genetics, 2005.
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
merge.analysis = function(pheno, pheno.col = 1, probs, K, addcovar, snps,
                 chr, start, end, model = c("additive", "dominance", "full"),
                 scan = c("one", "two"), output = c("lod", "bic"),
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
  pheno = pheno[!is.na(pheno[,pheno.col]),,drop = F]
  samples = intersect(rownames(pheno), dimnames(probs)[[1]])
  if(length(samples) == 0) {
    stop(paste("rownames(pheno) does not contain any sample IDs that",
         "match dimnames(probs)[[1]]. Please verify that both variables",
         "contain the same sample IDs."))
  } # if(length(samples) == 0)

  if(!missing(addcovar)) {
    addcovar = addcovar[rownames(addcovar) %in% samples,,drop = F]
    if(nrow(addcovar) == 0) {
      stop(paste("There are no common sample IDs between rownames(addcovar)",
           "and rownames(pheno). Please ensure that there is some overlap",
           "in sample IDs between rownames(addcovar) and rownames(pheno)."))
    } # if(nrow(addcovar) == 0)
    remove = which(rowSums(is.na(addcovar)) > 0)
    if(length(remove) > 0) {
      samples = samples[-remove]
      addcovar = addcovar[-remove,,drop = F]
    } # if(length(remove) > 0)
  } # if(!missing(addcovar))

  pheno = pheno[rownames(pheno) %in% samples,,drop = F]
  probs = probs[dimnames(probs)[[1]] %in% samples,,]
  probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,]
  K = K[rownames(K) %in% samples, colnames(K) %in% samples]
  K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), colnames(K))]
  print(paste("Mapping with", nrow(pheno), "samples."))

  stopifnot(all(rownames(pheno) == dimnames(probs)[[1]]))
  stopifnot(all(rownames(pheno) == rownames(K)))
  stopifnot(all(rownames(pheno) == colnames(K)))

  # Subset the SNPs to be the same.
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]

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
  sanger = lapply(sanger, strsplit, split = "\t")[[1]]
  sanger = matrix(unlist(sanger), length(sanger), length(sanger[[1]]),
           dimnames = list(sanger$POS, hdr), byrow = T)
  
  # Keep only high quality reads.
  qual.columns = grep("quality", colnames(sanger))
  sanger = sanger[apply(sanger[,qual.columns] == "1", 1, all),,drop = F]

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
  remove = na.rows[rowSums(sanger[na.rows,], na.rm = T) == 0]
  if(length(remove) > 0) {
    sanger = sanger[-remove,]
    ref = ref[-remove]
    pos = pos[-remove]
  } # if(length(remove) > 0)

  # Reorder the strains to match DOQTL founders.
  m = match(c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ",
              "PWK/PhJ", "WSB/EiJ"), colnames(sanger))
  if(any(is.na(m))) {
    stop("One of the dO founders is not in the SNP file.")
  } # if(any(is.na(m)))
  sanger = sanger[,m]
  rownames(sanger) = pos

  # Get the SDPs.
  sdps = data.frame(pos = as.numeric(pos),
         sdps = apply(sanger, 1, paste, collapse = ""), sanger)

  print("Calculating LOD")

  retval = NULL
  if(scan == "one") {
    retval = merge.scan1(pheno = pheno, pheno.col = pheno.col, probs = probs, 
             K = K, addcovar = addcovar, sanger = sanger, pos = pos,
             sdps = sdps, snps = snps, model = model, output = output)
  } else {
# TBD: implement pair scan.
    retval = merge.scan2(pheno, pheno.col, addcovar, sanger, pos, sdps, snps)
  } # else

  return(retval)

} # merge.analysis()


###
# Helper function for one-way scan.
# Using the BIC.
merge.scan1 = function(pheno, pheno.col, probs, K, addcovar, sanger, pos, sdps, 
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
    BIC.full = -2 * qtl.full$lod$lrs + (ncol(addcovar) +  dim(probs)[[2]]) * 
               log(nrow(pheno))
  } # if(output == "bic")
  
  stat = rep(0, nrow(sdps))

  # Loop through each MUGA SNP in the interval.
  index = 1
  for(i in 1:(nrow(snps)-1)) {

    # Get the Sanger SDPs between this MUGA SNP and the next.
    sdp.sub = sdps[sdps$pos >= snps[i,3] * 1e6 & sdps$pos < snps[i+1,3] * 1e6,]

    if(nrow(sdp.sub) > 0) {

      # Keep the unique SDPs.
      unique.sdp = sdp.sub[!duplicated(sdp.sub[,2]),]
      # Get the locations where the results should be placed.
      m = match(sdp.sub$sdps, unique.sdp$sdps)
      # Map the founder SNPs onto the DO genotypes.
      allele.probs = (0.5 * (probs[,,i] + probs[,,i+1])) %*% 
                      t(as.matrix(unique.sdp[,-1:-2]))

      lrs = matrixeqtl.snps(pheno = pheno[,pheno.col], geno = allele.probs, 
            K = err.cov, addcovar = addcovar)

      if(output == "bic") {
        BIC.reduced = -2 * lrs + (ncol(addcovar) + 3) * log(nrow(pheno))
        stat[index:(index + length(m) - 1)] = BIC.reduced[m]
      } else {
        stat[index:(index + length(m) - 1)] = (lrs / (2 * log(10)))[m]
      } # else
      index = index + length(m)

    } # if(nrow(sdp.sub) > 0)
  } # for(i)

  if(output == "bic") {
    bfull = rep(0, nrow(sdps))
    for(i in 1:(nrow(snps) - 1)) {
      bfull[sdps[,1] >= snps[i,3] * 1e6 & sdps[,1] <= snps[i+1,3] * 1e6] = BIC.full[i]
    } # for(i)

    return(data.frame(chr = rep(snps[1,2], nrow(sdps)), sdps, BIC.full = bfull, 
           BIC.red = stat))
  } else {
    return(data.frame(chr = rep(snps[1,2], nrow(sdps)), sdps, LOD = stat))
  } # else
  
} # merge.scan1


###
# Helper function for two-way scan.
merge.scan2 = function(pheno, pheno.col, addcovar, sanger, pos, sdps, snps) {

  stop("merge.scan2 not implemented yet.")

} # merge.scan2


###
# Plot for one-way merge analysis scan.
# Arguments: results: data.frame from merge.analysis.
# 
# 
merge.plot = function(results, 
  mgi.file = "ftp://ftp.jax.org/SNPtools/genes/MGI.20130703.sorted.txt.gz",
  highlight, highlight.col = "red", thr) {

  old.par = par(no.readonly = TRUE)

  start = results$pos[1] * 1e-6
  end   = results$pos[nrow(results)] * 1e-6

  bic = T
  if(colnames(results)[ncol(results)] == "LOD") {
    diff = results$LOD
    bic = F
  } else {
    diff = results$BIC.full - results$BIC.red
  } # else

  # Color the points with statistic > thr in red.
  col = rep(1, nrow(results))
  if(!missing(thr)) {
    if(bic) {
      points.ge.thr = which(results$BIC.full - results$BIC.red >= thr)
    } else {
      points.ge.thr = which(results$LOD >= thr)      
    } # else
    col[points.ge.thr] = 2
  } # else

  layout(matrix(1:2, 2, 1), heights = c(0.4, 0.6))
  par(plt = c(0.12, 0.99, 0, 0.9), las = 1)
  plot(results$pos * 1e-6, diff, ann = F, xaxt = "n", pch = 20, col = col)

  usr = par("usr")
  par(las = 3)
  if(bic) {
    mtext(side = 2, line = 3, text = "BIC.full - BIC.reduced")
  } else {
    mtext(side = 2, line = 3, text = "LOD of reduced model")
  } # else
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

} # merge.plot()


# Plot specific SDP locations along a chromosome.
# This idea is from Yalcin et.al., Genetics, 2005.
# results: the data.frame output from merge.analysis().
# sdps: strings of 0s and 1s that match some of the SDPs in column
#       2 of the merge.analysis() results.
# ...: additional arguments to pass to plot.
sdp.plot = function(results, sdps, ...) {

  plot(results[,1], rep(1, nrow(results)), col = 0, ylim = c(0, length(sdps)),
       yaxs = "i", ann = F, axes = F, ...)
  abline(h = 0:length(sdps))
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  mtext(side = 2, at = 1:length(sdps) - 0.5, text = sdps, line = 0.1, 
        adj = 1, las = 2)
  ax = axis(side = 1, labels = F)
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

