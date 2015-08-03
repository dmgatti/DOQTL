################################################################################
# Genome wide association mapping.
# We use parallel processing with multiple cores per node. We don't currently
# support multiple nodes.
# Daniel Gatti
# Dan.Gatti@jax.org
# Jan. 14, 2015
# Arguments: pheno: data.frame containing phenotypes in columns and samples in
#                   rows. Rownames must contain sample IDs.
#            pheno.col: the phenotype column to map.
#            probs: 3D numeric array containing the haplotype probabilities
#                   for each sample. Samples in rows, founders in columns,
#                   markers in slices. Samnple IDs must be in rownames. Marker
#                   names must be in dimnames(probs)[[3]].
#            K: List of kinship matrices, one per chromosome in markers.
#            addcovar: data.frame of additive covariates to use in the mapping.
#                      Samnple IDs must be in rownames.
#            markers: data.frame containing at least 3 columns with marker names
#                     chr, Mb postion.
#            sdp.file: character string containing the full path to the Sanger
#                      SDP file containing genomic positions and SDPs. This file
#                      is created using condense.sanger.snps().
#            ncl: The number of cores available.
################################################################################
# Contains scanone.assoc, s1.assoc, plot.scanone.assoc
################################################################################
scanone.assoc = function(pheno, pheno.col, probs, K, addcovar, markers,
                cross = c("DO", "CC", "HS"), sdp.file, ncl) {

  cl = makeCluster(ncl)
  registerDoParallel(cl)

  # Synch up markers and haplotype probs.
  markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]

  # Put the marker positions on a Mb scale.
  if(any(markers[,3] > 200)) {
    markers[,3] = markers[,3] * 1e-6
  } # if(any(markers[,3] > 200)

  # Synch up the sample IDs.
  samples = intersect(intersect(rownames(pheno), rownames(addcovar)),
            rownames(probs))
  pheno = pheno[samples,,drop = FALSE]
  addcovar = as.matrix(addcovar[samples,,drop = FALSE])
  probs = probs[samples,,]

  # Get the unique chromosomes.
  chr = unique(markers[,2])
  chr = lapply(chr, "==", markers[,2])
  chr = lapply(chr, which)
  names(chr) = unique(markers[,2])

  # Split up the data and create a list with elements for each 
  # chromosome.
  data = vector("list", length(chr))
  for(i in 1:length(chr)) {
    data[[i]] = list(pheno = pheno, pheno.col = pheno.col, 
                probs = probs[,,chr[[i]]], K = K[[i]], 
                addcovar = addcovar, markers = markers[chr[[i]],])
  } # for(i)
  names(data) = names(chr)

  rm(samples, probs, markers, K, addcovar)

  # Load the required libraries on the cores.
  clusterEvalQ(cl = cl, expr = library(DOQTL))
  clusterEvalQ(cl = cl, expr = library(Rsamtools))
  clusterEvalQ(cl = cl, expr = library(regress))

  res = foreach(obj = iter(data)) %dopar% {
    
    s1.assoc(obj, sdp.file)

  } # foreach(c)

  names(res) = names(data)

  stopCluster(cl)

  class(res) = c("scanone.assoc", class(res))
  return(res)

} # scanone.assoc()


# Helper function to perform association mapping on one autosome.
# obj: data object of the type created in scanone.assoc().
# sdp.file: Tabix file created by condense.sanger.snps().
s1.assoc = function(obj, sdp.file) {

  # Get the samples that are not NA for the current phenotype.
  sample.keep = which(!is.na(obj$pheno[,obj$pheno.col]))

  # Calculate variance components and change each kinship matrix to be the
  # error covariance correction matrix.
  mod = regress(obj$pheno[,obj$pheno.col] ~ obj$addcovar, ~obj$K,
        pos = c(TRUE, TRUE))
  obj$K = mod$sigma[1] * obj$K + mod$sigma[2] * diag(nrow(obj$K))
  rm(mod)

  # Read in the unique SDPs.
  tf = TabixFile(sdp.file)
  sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = obj$markers[1,2],
         ranges = IRanges(start = 0, end = 200e6)))[[1]]
  sdps = strsplit(sdps, split = "\t")
  sdps = matrix(unlist(sdps), ncol = 3, byrow = T)
  chr  = sdps[1,1]
  pos  = as.numeric(sdps[,2])
  sdps = as.numeric(sdps[,3])

  # Create a matrix of SDPs.
  sdp.mat = matrix(as.numeric(intToBits(1:2^8)), nrow = 32)
  sdp.mat = sdp.mat[8:1,]
  dimnames(sdp.mat) = list(LETTERS[1:8], 1:2^8)

  # Between each pair of markers, get the unique SDPs and their genomic
  # positions. Use the DO genoprobs to create DO genotypes.

  # Get a set of overlaps between the markers and SDP positions.
  sdp.gr = GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
  markers.gr = GRanges(seqnames = obj$markers[1,2], ranges = IRanges(
               start = c(0, obj$markers[,3]) * 1e6, 
               end = c(obj$markers[,3], 200) * 1e6))
  ol = findOverlaps(query = markers.gr, subject = sdp.gr)
  ol = split(subjectHits(ol), queryHits(ol))
  probs.idx = as.numeric(names(ol))
  unique.sdps = lapply(ol, function(z) { unique(sdps[z]) })
  num.sdps = sum(sapply(unique.sdps, length))
  geno = matrix(0, nrow = nrow(obj$pheno), ncol = num.sdps,
         dimnames = list(rownames(obj$pheno), 1:num.sdps))
  # This maps the positions in the geno matrix back to the genomic positions
  # of the SDPs.
  map = rep(0, length(pos))

  idx = 0
  # Start of chromosome, sdps before the first marker.
  i = 1
  rng = (idx + 1):(idx + length(unique.sdps[[i]]))
  geno[,rng] = obj$probs[,,probs.idx[i]] %*% sdp.mat[,unique.sdps[[i]]]
  map[ol[[i]]] = match(sdps[ol[[i]]], unique.sdps[[i]]) + idx
  idx = idx + length(unique.sdps[[i]])

  lrs = 0

  # SDPs bracketed by two markers.
  for(i in 2:(length(ol)-1)) {

    rng = (idx + 1):(idx + length(unique.sdps[[i]]))
    geno[,rng] = 0.5 * (obj$probs[,,probs.idx[i]] + obj$probs[,,probs.idx[i]+1]) %*% sdp.mat[,unique.sdps[[i]]]
    map[ol[[i]]] = match(sdps[ol[[i]]], unique.sdps[[i]]) + idx
    idx = idx + length(unique.sdps[[i]])

  } # for(i)

  # End of chromosome, sdps after the last marker.
  i = length(ol)
  rng = (idx + 1):(idx + length(unique.sdps[[i]]))
  geno[,rng] = obj$probs[,,probs.idx[i]] %*% sdp.mat[,unique.sdps[[i]]]
  map[ol[[i]]] = match(sdps[ol[[i]]], unique.sdps[[i]]) + idx

  # X Chromosome, separate females and males and then combine.
  if(chr == "X") {

    warning("X CHROMOSOME NOT FULLY IMPLEMENTED YET!")

    # Verify that sex is one of the covariates.
    if(!any("addcovar" == names(obj))) {

      stop(paste("In order to map on the X chromosome, you must supply the sex",
           "of each sample in \'addcovar\', even if they are all the same sex."))

    } # if(!any("addcovar" == names(obj)))

    if(length(grep("sex", colnames(obj$addcovar), ignore.case = T)) == 0) {

      stop(paste("In order to map on the X chromosome, you must supply the sex",
           "of each sample in \'addcovar\', even if they are all the same sex."))

    } # if(grep("sex", colnames(obj$addcovar), ignore.case = T))

    # Get the sex of each sample.
    sex.col = grep("sex", colnames(obj$addcovar), ignore.case = TRUE)
    females = which(obj$addcovar[,sex.col] == 0)
    males   = which(obj$addcovar[,sex.col] == 1)

    if(length(females) == 0 & length(males) == 0) {
      stop(paste("Sex is not coded using 0 for female and 1 for males. Please",
           "set the sex column in addcovar to 0 for females and 1 for males."))
    } # if(length(females) == 0 & length(males) == 0)

    # Fit a model with only sex.
    sst = sum((obj$pheno[sample.keep, obj$pheno.col] - mean(obj$pheno[sample.keep, obj$pheno.col]))^2)
    sex.rsq = lsfit(obj$addcovar[sample.keep,], obj$pheno[sample.keep, obj$pheno.col])
    sex.rsq = 1.0 - sum(sex.rsq$residuals^2) / sst

    # Separate the male and female genotypes in the same matrix. This will be 
    # a block matrix with all zeros for females in rows with male samples and
    # all zeros for males in rows with female samples.
    newgeno = matrix(0, nrow(geno), 2 * ncol(geno))
    newgeno[females,1:ncol(geno)] = geno[females,]
    newgeno[males,(ncol(geno)+1):ncol(newgeno)] = geno[males,]

    r2 = matrixeqtl.snps(pheno = obj$pheno[sample.keep,obj$pheno.col,drop = FALSE],
         geno = newgeno[sample.keep,,drop = FALSE],
         K = obj$K[sample.keep, sample.keep,drop = FALSE],
         addcovar = obj$addcovar[sample.keep,,drop = FALSE])
    lrs = r2[1:ncol(geno)] + r2[(ncol(geno)+1):length(r2)] - sex.rsq

  } else {

    # Autosomes.
    # Calculate the R^2 for each SDP.
    lrs = matrixeqtl.snps(pheno = obj$pheno[sample.keep,obj$pheno.col,drop = FALSE],
          geno = geno[sample.keep,,drop = FALSE],
          K = obj$K[sample.keep,sample.keep,drop = FALSE],
          addcovar = obj$addcovar[sample.keep,,drop = FALSE])

  } # else

  # Convert R^2 to LRS.
  lrs = -length(sample.keep) * log(1.0 - lrs)

  # Convert the LRS to p-values.
  pv = pchisq(2 * lrs, df = 1, lower.tail = FALSE)

  # Place the results in the correct locations and return.
  return(GRanges(seqnames = Rle(chr, length(pos)), ranges = IRanges(start = pos, 
         width = 1), p.value = pv[map]))

} # s1.assoc()


# Plotting for scanone.assoc.
plot.scanone.assoc = function(x, chr, bin.size = 10000, ...) {

  if(!missing(chr)) {
    x = x[names(x) %in% chr]
  } # if(!missing(chr))

  chrlen = get.chr.lengths()
  chrlen = chrlen[names(chrlen) %in% names(x)]
  chrlen = cumsum(chrlen)
  chrmid = c(0, chrlen[-length(chrlen)]) + diff(c(0, chrlen)) * 0.5
  names(chrmid) = names(chrlen)

  if(missing(chr)) {
    chr = names(x)
  } # if(!missing(chr))

  pos = lapply(x, start)
  pv =  lapply(x, mcols)
  pv =  lapply(pv, "[[", 1)

  for(c in 1:length(x)) {

    bins = seq(1, length(pv[[c]]), bin.size)
    bins[length(bins) + 1] = length(pv[[c]])
    pos2 = rep(0, length(bins))
    pv2  = rep(0, length(bins))

    for(i in 1:(length(bins)-1)) {
      wh = which.min(pv[[c]][bins[i]:bins[i+1]])
      wh = wh + bins[i] - 1
      pos2[i] = pos[[c]][wh] * 1e-6 + max(0, chrlen[c - 1])
      pv2[i]  = pv[[c]][wh]
    } # for(i)

    pos[[c]] = pos2
    pv[[c]]  = -log(pv2, 10)

  } # for(c)

  # If we are plotting more than one chormosome, color alternate 
  # chromosomes grey and black.
  col = 1
  chr = rep(names(x), sapply(pos, length))
  if(length(pos) > 1) {
    col = as.numeric(chr) %% 2 + 1
  } # if(length(chr) > 1)

  plot(unlist(pos), unlist(pv), pch = 16, col = c("black", "grey60")[col],
       las = 1, xaxt = "n", xlab = "", ylab = "-log10(p-value)", xaxs = "i", ...)
  mtext(text = names(chrmid), side = 1, line = 2.5, at = chrmid, cex = 2)

  if(length(pos) == 1) {

    axis(side = 1)

  } # if(length(pos) == 1)

} # plot.scanone.assoc()

