################################################################################
# Given a set of SNP locations, retrieve the full information, including
# strain calls and consequences.
# Daniel Gatti
# dan.gatti@jax.org
# May 26, 2015
################################################################################
# Arguments: gr: GRanges containing the SNP locations.
#            strains: character vector containing the strains to return.
#            snp.file: character vector containing the full path to the SNP
#                      file.
get.snp.details = function(gr, strains = character(), snp.file = 
                  "ftp://ftp.jax.org/sanger/current_snps/mgp.v4.snps.dbSNP.vcf.gz") {

  hdr = scanVcfHeader(snp.file)

  if(length(strains) == 0) {
    strains = samples(hdr)
  } # if(length(strains) == 0)

  param = ScanVcfParam(info = "CSQ", geno = c("GT", "FI"), 
          samples = strains, which = gr)

  snps = readVcf(file = snp.file, genome = "mm10", param = param)
  rowRanges(snps)$paramRangeID = names(rowRanges(snps))
  snps = genotypeCodesToNucleotides(snps)
  snps = expand(snps)

  alt = sapply(fixed(snps)$ALT, as.character)

  # Try to convert the consequences into one line.
  csq = info(snps)$CSQ
  keep = rep(FALSE, nrow(snps))
  for(i in 1:nrow(snps)) {

    spl = strsplit(csq[[i]], split = "\\|")
    spl = matrix(unlist(spl), ncol = length(spl))
    spl = spl[,spl[1,] == alt[i], drop = FALSE]
    
    tmp = paste0(spl[2,], spl[5,])
    spl = spl[,!duplicated(tmp), drop = FALSE]
    csq[[i]] = apply(spl, 2, paste, collapse = "|")

    keep[i] = alt[i] %in% unique(substring(sub("/", "", geno(snps)$GT[i,]), 1, 1))

  } # for(i)

  snps = snps[keep]

  snps

} # get.snp.details()



# Test code
#library(GenomicRanges)
#setwd("D:/234_DO_Harrison/QTL/assoc_loco_lm/hq_snps/")
#load("lifespan_assoc_loco_lm_hq_chr18_zoom.Rdata")

#assoc = assoc[-log10(assoc[,12]) > 12,]
#gr = GRanges(seqnames = assoc[,1], range = IRanges(start = assoc[,2],
#     width = 1), mcols = assoc[,-(1:2)])
#load("C:/Users/dgatti/Documents/ensembl/Mus_musculus.GRCm38.77.Rdata")
#gr = subsetByOverlaps(gr, ensembl)

#strains = c("A_J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")

#x = get.snp.details(gr, strains)
