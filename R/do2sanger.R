################################################################################ 
# Read in the genotypes of DO mice, assign the most probable genotype based on 
# the posterior probabilities and map the Sanger SNPs onto the DO samples.
# These functions split the DO genotypes at the midpoint between MUGA markers.
# Daniel Gatti
# Dan.Gatti@jax.org
# Sept. 7, 2013
################################################################################
# Arguments: do.files: character vector of *.genotype.probs.Rdata files that 
#                      contain the posterior probabilities.
#            snps: data.frame containing 4 columns containing the SNP IDs,
#                  the chromosome, the Mb location and the cM location.
#            snp.file: character containing the full path to the Sanger SNP file.
#                      Defaults to the SNP files at ftp.jax.org.
do2sanger = function(do.files, snps, output.file = "do2sanger.txt", 
            snp.file = "ftp://ftp.jax.org/SNPtools/variants/cc.snps.NCBI38.txt.gz") {
  if(missing(do.files)) {
    stop("do.files cannot be missing. Please enter a set of DO sample files.")
  } # if(missing(do.files))
  if(missing(snps)) {
    stop(paste("snps cannot be missing because we need the locations of the",
         "SNPs. Please provide SNP locations that match the SNPS in the",
         "do.files."))
  } # if(missing(do.files))
  # Read in the first file to get the genotype and SNP dimensions.
  load(do.files[1])
  states = colnames(prsmth)
  snps = snps[snps[,1] %in% rownames(prsmth),]
  prsmth = prsmth[rownames(prsmth)[[3]] %in% snps[,1],]
  gt = matrix("", length(do.files), nrow(prsmth), dimnames = list(
       sub("\\.genotype\\.probs\\.Rdata$", "", do.files), rownames(prsmth)))
  # Read in the DO files and call the maximum genotypes.
  print("Calling DO genotypes...")
  for(i in 1:length(do.files)) {
    load(do.files[i])
    gt[i,] = states[apply(prsmth, 1, which.max)]
  } # for(i)
  # Split up the MUGA markers by chromosomes.
  snps = snps[snps[,1] %in% rownames(prsmth),]
  snps = split(x = snps, f = snps[,2])
  sanger = get.variants(file = snp.file, chr = 1, start = 0, end = 4,
           strains = c("129S1/SvImJ", "A/J", "C57BL/6J", "CAST/EiJ", 
           "NOD/ShiLtJ", "NZO/HlLtJ", "PWK/PhJ", "WSB/EiJ"))
  keep = setdiff(1:ncol(sanger), grep("quality", colnames(sanger)))
  for(c in 1:length(snps)) {
    curr.chr = names(snps)[c]
    print(paste("CHR", curr.chr))
    num.snps = 0
    # Start writing SNPs at the beginning of the chromosome 1.
    print("Mapping Sanger SNPs to DO samples...")
    con = file(description = paste("chr", curr.chr, ".", output.file, sep = ""),
          open = "w")
    writeLines(paste("ID", "CHR", "GRC38.bp", "Ref", "Alt", 
               paste(rownames(gt), collapse = "\t"), sep = "\t"), con)
    gc()
    # Start of Chr.
    result = do2sanger.helper(snp.file = snp.file, chr = curr.chr, start = 0, 
             end = mean(snps[[curr.chr]][1:2,3]), geno = gt[,snps[[curr.chr]][1,1]],
             keep = keep)
    if(!is.null(result)) {
      result = apply(result, 1, paste, collapse = "\t")
      writeLines(result, con)
      num.snps = num.snps + nrow(result)
    } # if(!is.null(result))
    # Middle SNPs.
    for(i in 2:(nrow(snps[[curr.chr]]) - 1)) {
      if(i %% 100 == 0) {
        print(i)
      } # if(i %% 10 == 0)
      if(all(!is.na(snps[[curr.chr]][(i-1):(i+1),3]))) {
        result = do2sanger.helper(snp.file = snp.file, chr = curr.chr,
                 start = mean(snps[[curr.chr]][(i-1):i,3]), 
                 end = mean(snps[[curr.chr]][i:(i+1),3]), 
                 geno = gt[,snps[[curr.chr]][i,1]], keep = keep) 
        if(!is.null(result)) {
          result = apply(result, 1, paste, collapse = "\t")
          writeLines(result, con)      
          num.snps = num.snps + nrow(result)
        } # if(!is.null(result))
      } # if(all(!is.na(snps[[curr.chr]][(i-1):(i+1),3])))
    gc()
    } # for(i)
    # End of Chr.
    if(all(!is.na(snps[[curr.chr]][(i-1):(i+1),3]))) {
      result = do2sanger.helper(snp.file = snp.file, chr = curr.chr,
               start = mean(snps[[curr.chr]][(nrow(snps[[curr.chr]])-1):nrow(snps[[curr.chr]]),3]), 
               end = 200, geno = gt[,snps[[curr.chr]][nrow(snps[[curr.chr]]),1]], keep = keep)
      if(!is.null(result)) {
        result = apply(result, 1, paste, collapse = "\t")
        writeLines(result, con)
        num.snps = num.snps + nrow(result)
      } # if(!is.null(result))
    } # if(all(!is.na(snps[[curr.chr]][(i-1):(i+1),3])))
    close(con)
    gc()
    print(paste("CHR", curr.chr, ":", num.snps))
  } # for(c)
} # do2sanger()
##########
# Helper function to get Sanger SNPs and place them in the DO samples.
do2sanger.helper = function(snp.file, chr, start, end, geno, keep) {
  sanger = get.variants(file = snp.file, chr = chr, start = start, end = end,
           strains = c("129S1/SvImJ", "A/J", "C57BL/6J", "CAST/EiJ", 
           "NOD/ShiLtJ", "NZO/HlLtJ", "PWK/PhJ", "WSB/EiJ"), quality = 1)
  if(length(sanger) > 0 && nrow(sanger) > 0) {
    sanger = sanger[,keep]
    geno.cols = (ncol(sanger) - 7):ncol(sanger)
    map = do.colors[,1:2]
    map = map[match(colnames(sanger)[geno.cols], map[,2]),]
    colnames(sanger)[geno.cols] = map[,1]
  
    # Here, we assume that the founders are homozygous. This may not be the case 
    # and this is an area for improvement. 
    sanger[,geno.cols] = apply(sanger[,geno.cols], 2, substring, 1, 1)
    hdr = sanger[,1:5]
    sanger = as.matrix(sanger[,-1:-5])
    # Turn the sanger SNPs into 0 for reference allele and 1 for alternate allele.
    sanger = matrix(as.numeric(sanger != hdr[,4]), nrow(sanger), ncol(sanger),
             dimnames = dimnames(sanger))
    # Create a 2 x n matrix with founder haplotypes for each DO founder.
    samples = names(geno)
    geno = matrix(unlist(strsplit(geno, split = "")), nrow = 2, dimnames = 
           list(NULL, samples))
    # Convert the founder haplotypes to numeric genotypes.
##### NOTE: We are assuming bi-allelic SNPs here! This may not be the case
#####       for all Sanger SNPs. 
    geno = matrix(sanger[,geno[1,]] + sanger[,geno[2,]],
           nrow(sanger), ncol(geno), dimnames = list(NULL, samples))
    # Code the major allele as 0 .
    num.zero = rowMeans(geno == 0, na.rm = TRUE)
    num.two  = rowMeans(geno == 2, na.rm = TRUE)
    swap.rows = which(num.two > num.zero)
    geno[swap.rows,] = 2 - geno[swap.rows,]
    return(cbind(hdr, geno)) 
  } else {
    return(NULL)
  } # else
} # do2sanger.helper()
