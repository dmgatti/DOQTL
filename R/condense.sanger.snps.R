################################################################################
# Get the Sanger SNPs and condense them down to unique strain distribution
# patterns at each location.
# This will write out a tab delimited file with the chr, position and SDP
# code. The code will be binary with the founders listed in alphabetical
# order, i.e. A-H. A '1' indicates that the strain contains an alternate allele
# and '0' indicates the reference. The binary code run from right to left.
# i.e. '19' is 00010011.
################################################################################
### NOTE: we consider heterozygotes and trimorphic SNPs as bimorphic. ##########
################################################################################
# Daniel Gatti
# Dan.Gatti@jax.org
# Dec. 6, 2014
################################################################################
# Create a function to generate the file, given the array and SNP file.
condense.sanger.snps = function(markers, snp.file, strains, out.file, ncl) {

  # Get the header information.
  hdr = scanVcfHeader(snp.file)

  # Check that the strains are in the file.
  if(!all(strains[strains != "C57BL_6J"] %in% samples(hdr))) {
    stop(paste("Some strains are not in the SNP file. Strains in file:",
         paste(samples(hdr), collapse = ",")))
  } # if(!all(strains %in% samples(hdr)))

  # Set up for parallel execution.
  cl = NULL
  if(!missing(ncl)) {

    cl = makeCluster(spec = ncl, type = "MPI", outfile = "cluster_out.txt")
    registerDoParallel(cl = cl)
    clusterEvalQ(cl = cl, expr = library(DOQTL))
    clusterEvalQ(cl = cl, expr = library(Rsamtools))
    clusterExport(cl = cl, "markers")

  } # if(!missing(ncl))

  # Get the SNPs and each chromosome.
  files = foreach(chr = iter(seqlevels(markers))) %dopar% {

    print(chr)

    # Get the markers on the current chromosome.
    cur.chr = markers[seqnames(markers) == chr]
    
    # Get SNPs from the beginning of the chromosome to the first marker.
    gr = GRanges(seqnames = seqnames(cur.chr)[1], ranges = IRanges(start = 0,
                 end = 200e6))

    sdps = get.snp.patterns(snp.file = snp.file, gr = gr,
           strains = strains, polymorphic = TRUE, filter.by.qual = TRUE)

    sdps = paste(sdps[,1], sdps[,2], sdps[,ncol(sdps)], sep = "\t")

    filename = paste0("chr", chr, ".txt")
    outf = file(description = filename, open = "w")
    writeLines(sdps, con = outf, sep = "\n")
    close(outf)
    filename

  } # foreach(c)

  if(!missing(ncl)) {

    stopCluster(cl)

  } # if(!missing(ncl))

  # Read in each of the files and write them out to one, large file.
  files = unlist(files)
  outf = file(description = out.file, open = "w")
  writeLines("# DO Sanger v4 SDPs, high quality, polymorphic SNPs only. DMG 4/3/2015", 
             con = outf, sep = "\n")
  writeLines(paste0("#", paste(strains, sep = "\t")), con = outf, sep = "\n")
  writeLines("#chr\tGRCm38\tSDP", con = outf, sep = "\n")

  # We have to change the scientific notation penalty to print out the large
  # integers without using scientific notation, which kills the Tabix indexing.
  old.op = options("scipen")
  options(scipen = 10)

  for(i in 1:length(files)) {

    x = read.delim(files[i], header = FALSE)
    x = paste(x[,1], x[,2], x[,3], sep = "\t")
    writeLines(x, con = outf, sep = "\n")

  } # for(i)

  close(outf)

  options(scipen = old.op)

  bgzip(file = out.file, overwrite = TRUE)
  indexTabix(file = paste0(out.file, ".bgz"), seq = 1, start = 2, end = 2)

} # condense.sanger.snps()


# Function to read the Sanger SNP VCF and produce SNP patterns.
# Filter to keep only PASS SNPs and SNPs with all '1' quality.
get.snp.patterns = function(snp.file = snp.file, gr = gr,
                   strains = strains, polymorphic = TRUE, filter.by.qual = TRUE) {

  # Get the header information.
  hdr = scanVcfHeader(snp.file)

  # Check that the strains are in the file.
  if(!all(strains[strains != "C57BL_6J"] %in% samples(hdr))) {
    stop(paste("Some strains are not in the SNP file. Strains in file:",
         paste(samples(hdr), collapse = ",")))
  } # if(!all(strains %in% samples(hdr)))

  # Check for C57BL/6J.
  bl6 = which(strains == "C57BL_6J")
  new.strains = strains
  if(length(bl6) > 0) {
    new.strains = strains[-bl6]
  } # if(any(is.na(keep)))

  # Get the SNPs on this chromosome.
  param = ScanVcfParam(geno = c("GT", "FI"), samples = strains, 
          which = gr)
  snps = readVcf(file = snp.file, genome = "mm10", param = param)

  message(paste("Retrieved", nrow(snps), "SNPs."))

  # If the user wants to filter by quality, keep only SNPs where all calls have
  # quality == 1.
  if(filter.by.qual) {

    snps = snps[rowSums(geno(snps)$FI, na.rm = TRUE) == ncol(snps)]
    message(paste("Retaining", nrow(snps), "high quality SNPs."))

  } # if(filter.by.qual)

  # If the user requested only polymorphic SNPs, remove non-polymorphic ones.
  if(polymorphic) {

    snps = snps[rowSums(geno(snps)$GT == "0/0", na.rm = TRUE) < ncol(snps)]
    message(paste("Retaining", nrow(snps), "polymorphic SNPs."))

  } # if(polymorphic)

  # Get numeric allele calls.
  calls = (geno(snps)$GT != "0/0") * 1

  # Add in C57BL/6J, if it was requested.
  if(length(bl6) == 1) {

    calls = cbind(calls[,1:(bl6-1),drop = FALSE], C57BL_6J = rep(0, nrow(calls)), 
                  calls[,(bl6):ncol(calls),drop = FALSE])
    stopifnot(colnames(calls) == strains)

  } # if(length(bl6) == 1)

  # Flip the alleles for SDPs that have more ones than zeros.
  # NOTE: We treat all SNPs as bimorphic!
  flip = which(rowSums(calls, na.rm = TRUE) > ncol(calls) / 2)
  calls[flip,] = 1 - calls[flip,]

  alt = CharacterList(alt(snps))
  alt = unstrsplit(alt, sep = ",")
  snps = data.frame(chr = as.vector(seqnames(snps)), pos = start(snps), 
         id = names(rowRanges(snps)), ref = as.vector(ref(snps)), alt = alt, calls, 
         sdp = calls %*% 2^((ncol(calls)-1):0), stringsAsFactors = FALSE)

  # Fix the 129S1 column name.
  colnames(snps) = sub("^X", "", colnames(snps))

  return(snps)

} # get.snp.patterns()

