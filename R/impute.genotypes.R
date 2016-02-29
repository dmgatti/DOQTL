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
# outfile: Character string containing the path and filename to write out to.
#          Should end with "Rdata".
impute.genotypes = function(gr, probs, vcf.file, hq = TRUE, 
                   cross = c("DO", "CC", "DOF1", "HS", "HSrat", "other"),
                   filename = "imputed_snps.Rdata") {

  if(dim(probs)[3] != nrow(markers)) {
    print(paste0("The dim(probs)[3] (", dim(probs)[3], "must equal nrow(markers) (",
          nrow(markers), ")."))
  } # if(dimnames(probs)[3] != nrow(markers))

  if(any(dimnames(probs)[[3]] != markers[,1])) {
    print(paste0("The dimnames(probs)[[3]] must equal nrow(markers)."))
  } # if(any(dimnames(probs)[[3]] != markers[,1]))

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
  wh = which(as.character(markers[,2]) == as.character(seqnames(gr)) & markers[,3] >= start(gr) & 
       markers[,3] <= end(gr))
  markers = markers[c(wh[1] - 1, wh, wh[length(wh)] + 1),]
  probs = probs[,,markers[,1]]

  brks = cut(pos, markers[,3])
  pos = paste(markers[1,2], pos, sep = "_")
  pos = split(pos, brks)
  nr = nrow(mat)
  brks2 = rep(brks, each = 8)
  mat = split(mat, brks2)
  rm(brks, brks2)
  mat = lapply(mat, matrix, nrow = nr)
  
  outfile = file(description = filename, open = "w")
  writeLines(text = paste("SNP", paste(rownames(probs), sep = " "), sep = " "),
             con = outfile)

  locs = sapply(mat, ncol)
  snps = matrix(0, nrow = sum(locs), ncol = nrow(probs), dimnames =
         list(1:sum(locs), rownames(probs)))
  locs = c(0, cumsum(locs))

  for(i in 1:(length(mat)-1)) {

    pr = apply(probs[,,i:(i+1)], 1:2, mean, na.rm = TRUE)
    s = round(2 * pr %*% mat[[i]])
    
    stopifnot(range(s) == c(0,2))
    
    rng = (locs[i] + 1):locs[i+1]
    snps[rng,] = t(s)
    rownames(snps)[rng] = pos[[i]]

  } # for(i)

  i = length(mat)
  pr = probs[,,i]
  s = round(2 * pr %*% mat[[i]])
    
  stopifnot(range(s) == c(0,2))
    
  rng = (locs[i] + 1):locs[i+1]
  snps[rng,] = t(s)
  rownames(snps)[rng] = pos[[i]]

  print(paste("Writing out file:", filename))
  outfile = file(description = filename, open = "w")
  saveRDS(snps, file = outfile)
  flush(con = outfile)
  close(con = outfile)
  
  print(showConnections(all = TRUE))

  # We need this or else R returns the file from the saveRDS() function
  # and keeps the connection open.
  return()
  
} # impute.genotypes()
