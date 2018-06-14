################################################################################
# Impute the Sanger genotypes onto the DO genomes using the haplotype 
# probabilities.
# Daniel Gatti
# dan.gatti@jax.org
# Feb. 18, 2016
################################################################################
# Arguments: 
# gr: GRanges object that contains the genomic range in which to impute SNPs.
#     For now, this must be a single, continuous range.
# probs: 3D numeric array containing the haplotype probabilities. samples x 
#        8 founders x markers. All dimensions must be contain dimnames.
# markers: data.frame containing at least 3 columns that include the marker
#          ID, chr and postion of each marker. nrow must be the same as
#          dim(probs)[3] and markers[,1] must equal dimnames(probs)[[3]].
# vcf.file: String containing the full path to the Sanger SNP VCF file.
# hq: Boolean indicating whether to use only high quality SNPs. Default = TRUE.
# cross: Character string that is the cross type. One of "DO", "CC", "DOF1", 
#        "HS", "HSrat", "other")
impute.genotypes = function(gr, probs, markers, vcf.file, hq = TRUE, 
                   cross = c("DO", "CC", "DOF1", "HS", "HSrat", "other")) {
  if(dim(probs)[3] != nrow(markers)) {
    print(paste0("The dim(probs)[3] (", dim(probs)[3], "must equal nrow(markers) (",
          nrow(markers), ")."))
  } # if(dimnames(probs)[3] != nrow(markers))
  if(any(dimnames(probs)[[3]] != markers[,1])) {
    print(paste0("The dimnames(probs)[[3]] must equal nrow(markers)."))
  } # if(any(dimnames(probs)[[3]] != markers[,1]))
  # Synch up the marker locations.
  probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]
  markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
  stopifnot(markers[,1] == dimnames(probs)[[3]])
  # Some of the MUGA series markers occur at the same position.
  # Keep only unique ones.
  unique.markers = which(!duplicated(markers[,3]))
  markers = markers[unique.markers,]
  probs = probs[,,unique.markers]
  rm(unique.markers)
  # Check the GRanges ranges to see if they're in bp or Mb. We need to change 
  # them to bp for scan VCF to work.
  if(max(start(gr), end(gr)) < 200) {
    gr = GRanges(seqnames = seqnames(gr), range = IRanges(
         start = start(gr) * 1e6, end = end(gr) * 1e6))
  } # if(max(start(gr), end(gr)) < 200)
  # Put the markers on a bp scale.
  if(max(markers[,3]) < 300) {
    markers[,3] = markers[,3] * 1e6
  } # if(max(markers[,3] < 300)
  # Get the VCF header information.
  hdr = scanVcfHeader(vcf.file)
  # Get the founder samples to extract.
  samples = NULL
  if(cross == "DO" | cross == "CC") {
     samples = c("A_J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ",
                 "PWK_PhJ", "WSB_EiJ")
  } else {
    stop("Crosses other than the CC or DO have not bee implemented yet.")
  } # else
  # Create a VCF parameters.
  param = ScanVcfParam(samples = samples, fixed = c("ALT", "FILTER", "QUAL"),
          geno = c("GT", "FI"), which = gr)
  sanger = readVcf(file = vcf.file, genome = "mm10", param = param)
  snps = NULL
  if(nrow(sanger) > 0) {
    # Keep only high quality SNPs if requested.
    if(hq) {
      sanger = sanger[rowRanges(sanger)$FILTER == "PASS",]
    } # if(hq)
    mat = genotypeToSnpMatrix(sanger)
    mat = matrix(as.numeric(mat$genotypes), nrow = nrow(mat$genotypes), 
          ncol = ncol(mat$genotypes), dimnames = dimnames(mat$genotypes))
    # Add C57BL/6J to the SNP matrix.
    if(cross == "DO" | cross == "CC") {
      mat = rbind(mat[1,,drop = FALSE], C57BL_6J = rep(1, ncol(mat)), mat[2:7,])
    } # if(cross = "DO" | cross == "CC")
    pos = start(sanger)
    rm(sanger)
    # Keep only polymorphic SNPs.
    keep = which(colSums(mat == 1) < nrow(mat))
    mat = mat[,keep]
    pos = pos[keep]
    # Keep only bimorphic SNPs.
    x = sapply(apply(mat, 2, unique), length)
    keep = which(x == 2)
    mat = mat[,keep]
    pos = pos[keep]
    rm(x)
    # Convert SNPs to 0s and 1s.
    mat = (mat != 1) * 1
    # Impute the SNPs.
    wh = which(as.character(markers[,2]) == as.character(seqnames(gr)))
    markers = markers[wh,]
    probs = probs[,,markers[,1]]
    wh = which(markers[,3] >= start(gr) & markers[,3] <= end(gr))
    wh = unique(c(max(1, wh[1] - 1), wh, min(wh[length(wh)] + 1, nrow(markers))))
    markers = markers[wh,]
    probs = probs[,,markers[,1]]
    # Make breakpoints between markers and get the unique SDPs between each 
    # pair of markers.
    brks = cut(pos, markers[,3])
    pos = paste(markers[1,2], pos, sep = "_")
    pos = split(pos, brks)
    nr = nrow(mat)
    brks2 = rep(brks, each = 8)
    mat = split(mat, brks2)
    rm(brks, brks2)
    mat = lapply(mat, matrix, nrow = nr)
    mat = lapply(mat, t)
    locs = sapply(mat, nrow)
    snps = matrix(0, nrow = sum(locs), ncol = nrow(probs), dimnames =
           list(1:sum(locs), rownames(probs)))
    locs = c(0, cumsum(locs))
    mat.gt.0 = which(sapply(mat, length) > 0)
    for(i in mat.gt.0[1:(length(mat.gt.0)-1)]) {
      print(i)
      # The range of rows in snps to populate.
      rng = (locs[i] + 1):locs[i+1]
      # We sum the probs at the two surrounding markers,
      # but we DON'T divide by 2 because we want homozygotes
      # to be 0 or 2 and hets to be 1.
      pr = probs[,,i] + probs[,,(i+1)]
      snps[rng,] = round(tcrossprod(mat[[i]], pr))
      rownames(snps)[rng] = pos[[i]]
    } # for(i)
    i = mat.gt.0[length(mat.gt.0)]
    rng = (locs[i] + 1):locs[i+1]
    pr = 2 * probs[,,i]
    snps[rng,] = round(tcrossprod(mat[[i]], pr))
    rownames(snps)[rng] = pos[[i]]
    stopifnot(range(snps) == c(0,2))
  } # if(nrow(sanger) > 0)
  return(snps)
  
} # impute.genotypes()
