# Association mapping unit tests.
test_retrieve_correct_snps = function() {

  library(VariantAnnotation)

  snp.file = "ftp://ftp.jax.org/sanger/current_snps/mgp.v4.snps.dbSNP.vcf.gz"

  gr = GRanges(seqnames = 6, ranges = IRanges(start = 122e6, end = 123e6))
  sanger = get.snp.patterns(snp.file = snp.file, gr = gr,
           strains = c("A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ",
           "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ"), 
           polymorphic = TRUE, filter.by.qual = TRUE)
  rownames(sanger) = 1:nrow(sanger)

  # Do the same by hand.
  hdr = scanVcfHeader(snp.file)
  strains = samples(hdr)[c(5, 2, 21, 23, 13, 25, 28 )]
  param = ScanVcfParam(geno = c("GT", "FI"), samples = strains,
          which = gr)
  sanger2 = readVcf(file = snp.file, genome = "mm10", param = param)

  # Filter by quality.
  sanger2 = sanger2[rowSums(geno(sanger2)$FI, na.rm = T) == ncol(sanger2)]

  # Keep only polymorphic loci.
  sanger2 = sanger2[rowSums(geno(sanger2)$GT == "0/0") < ncol(sanger2)]

  # Get genotypes and convert to numeric values.
  geno = geno(sanger2)$GT
  geno = cbind(geno[,1], rep("0/0", nrow(geno)), geno[,2:7])
  geno = (geno != "0/0") * 1
  colnames(geno) = c("A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ",
                     "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")

  # Flip the bits for SDPs with more than 4 1s.
  flip = which(rowSums(geno, na.rm = T) > 4)
  geno[flip,] = 1 - geno[flip,]
  range(rowSums(geno, na.rm = T))

  sdps = geno %*% 2^(7:0)

  alt = CharacterList(alt(sanger2))
  alt = unstrsplit(alt, sep = ",")
  sanger2 = data.frame(CHROM = as.vector(seqnames(sanger2)), POS = start(sanger2),
            ID = names(rowRanges(sanger2)), REF = as.vector(ref(sanger2)),
            ALT = alt, geno, sdp = sdps, stringsAsFactors = FALSE)
  rownames(sanger2) = 1:nrow(sanger2)

  for(i in 1:ncol(sanger)) {
    checkEquals(target = sanger2[,i], current = sanger[,i], 
         msg = "get.snp.patterns didn't return the correct values.")
  } # for(i)

} # test_retrieve_correct_snps()


