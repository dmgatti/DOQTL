################################################################################
# Extract the intensity or genotype data from the raw GeneSeek data files.
# Arguments: in.path: character vector of full paths to the GeneSeek data
#                     directories.
#            prefix: character vector of same length as in.path containing
#                    a prefix to add to each sample ID in data sets being 
#                    processed.
#            out.path: character, full path to output directory.
# Output: Writes out call.rate.batch.txt and genotypes.txt for all options.
#         Writes out x.txt and y.txt when type = "intensity".
################################################################################
extract.raw.data = function(in.path = ".", prefix, out.path = ".", 
                   array = c("gigamuga", "megamuga", "muga")) {

  array = match.arg(array)

  if(!missing(prefix)) {
    if(length(in.path) != length(prefix)) {
      stop(paste("extract.raw.data: The 'in.path' and 'prefix' matrices must",
           "be the same length. Each prefix will be added to the matching",
           "data set in in.path."))
    } # if(length(in.path) != length(prefix))
  } # if(!missing(prefix))

  tmp = sub("/$", "", in.path)
  if(!all(file.exists(tmp))) {
    stop(paste("extract.raw.data: Some of the in.path directories do not exist.",
	     paste(in.path[!file.exists(tmp)], collapse = ",")))
  } # if(!all(file.exists(in.path)))
  
  snps = NULL
  if(array == "muga") {
    # We need the next line to satisfy R CMD check --as-cran
    muga_snps = NULL
    load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
    snps = muga_snps
  } else if(array == "megamuga") {
    # We need the next line to satisfy R CMD check --as-cran
    MM_snps = NULL
    load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
    snps = MM_snps
  } else if(array == "gigamuga") {
    # We need the next line to satisfy R CMD check --as-cran
    GM_snps = NULL
    load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
    snps = GM_snps
   } # else

  # Write out headers for the files.  This will overwrite existing files.
  x_file = file(paste(out.path, "x.txt", sep = "/"), open = "w")
  y_file = file(paste(out.path, "y.txt", sep = "/"), open = "w")
  g_file = file(paste(out.path, "geno.txt", sep = "/"), open = "w")
  writeLines(text = snps[-length(snps[,1]),1], con = x_file, sep = "\t")
  writeLines(text = snps[length(snps[,1]),1],  con = x_file, sep = "\n")
  writeLines(text = snps[-length(snps[,1]),1], con = y_file, sep = "\t")
  writeLines(text = snps[length(snps[,1]),1],  con = y_file, sep = "\n")
  writeLines(text = snps[-length(snps[,1]),1], con = g_file, sep = "\t")
  writeLines(text = snps[length(snps[,1]),1],  con = g_file, sep = "\n")

  call.rate.batch = NULL
  for(i in 1:length(in.path)) {

    # Get the sample IDs from the Sample_Map.txt file.
    samplefile = dir(path = in.path[i], pattern = "Sample_Map.txt", full.names = TRUE)

    # If not found, then quit.
    if(length(samplefile) == 0) {
      stop(paste("No file called 'Sample_Map.txt' was found in directory",
           in.path[i], ".  Please make sure that the Sample_Map file is unzipped and",
           "in the specified directory."))
    } # if(length(samplefile) == 0)

    samples = read.delim(samplefile, stringsAsFactors = FALSE)$Name
    samples = samples[nchar(samples) > 0]

    # Find a file with "FinalReport" in the filename.
    rawfile = dir(path = in.path[i], pattern = "FinalReport", full.names = TRUE)
    rawfile = rawfile[grep("txt", rawfile)]

    # If not found, then quit.
    if(length(rawfile) == 0) {
      stop(paste("No file with 'FinalReport' in the filename was found in directory",
           in.path[i], ".  Please make sure that the FinalReport file is unzipped and",
           "in the specified directory."))
    } # if(length(rawfile) == 0)

    # If there is more than one FinalReport file, then quit.
    if(length(rawfile) > 1) {
      stop(paste("There is more than one file with FinalReport in the filename.",
           "Please place only one data set in each directory."))
    } # if(length(rawfile) > 1)

    # Read in the first sample.  The current format requires us to skip 9 lines
    # because there is no comment delimiter at the top of the file.
    print(paste("Reading", rawfile, "..."))
    rawfile = file(rawfile, open = "r")
    data = readLines(con = rawfile, n = 10)
    hdr = strsplit(data, split = "\t")
    cn = hdr[[10]]

    # Verify that we have all of the column names that we expect.
    column.names = c("SNP Name", "Sample ID", "X", "Y", "Allele1 - Forward",
                     "Allele2 - Forward")
    columns = match(column.names, cn)  
    if(any(is.na(columns))) {
      stop(paste("All of the expected column names were not found in the",
           "FinalReport file. The missing column(s) are:", paste(
           colnames[is.na(columns)], collapse = ",")))
    } # if(any(is.na(columns)))
    cr = rep(0, length(samples))
    samples.in.data = rep("", length(samples))

    # We read the files in and write them out to conserve memory. As computers
    # get larger, we may be able to keep everything in memory.
    for(j in 1:length(samples)) {
      # Read in the data for one sample.
      data = readLines(con = rawfile, n = nrow(snps))
      data = strsplit(data, split = "\t")
      data = matrix(unlist(data), nrow = length(data[[1]]), ncol = nrow(snps))
      dimnames(data) = list(cn, data[1,])
      samples.in.data[j] = data[rownames(data) == "Sample ID",1]
      print(paste("Sample", j, "of", length(samples), ":", samples.in.data[j]))
      if(!missing(prefix)) {
        samples.in.data[j] = paste(prefix[i], samples.in.data[j], sep = "")
      } # if(!missing(prefix)) 
      # Sort the data to match the SNP order.
      data = data[,match(snps[,1], colnames(data))]
      # X
      writeLines(samples.in.data[j], con = x_file, sep = "\t")
      xint = data[rownames(data) == "X",]
      writeLines(xint[-length(xint)], con = x_file, sep = "\t")
      writeLines(xint[length(xint)], con = x_file)
      # Y
      writeLines(samples.in.data[j], con = y_file, sep = "\t")
      yint = data[rownames(data) == "Y",]
      writeLines(yint[-length(yint)], con = y_file, sep = "\t")
      writeLines(yint[length(yint)], con = y_file)
      # Genotype
      geno = paste(data[rownames(data) == "Allele1 - Forward",],
             data[rownames(data) == "Allele2 - Forward",], sep = "")
      cr[j] = mean(geno != "--")
      alleles = unique(geno)
      hets = alleles[!alleles %in% c("--", "AA", "CC", "GG", "TT")]
      geno[geno == "AA"] = "A"
      geno[geno == "CC"] = "C"
      geno[geno == "GG"] = "G"
      geno[geno == "TT"] = "T"
      geno[geno == "--"] = "-"
      geno[geno %in% hets] = "H"
 
      writeLines(samples.in.data[j], con = g_file, sep = "\t")
      writeLines(geno[-length(geno)], con = g_file, sep = "\t")
      writeLines(geno[length(geno)], con = g_file)
    } # for(j)
    close(rawfile)
    names(cr) = samples.in.data
    call.rate.batch = rbind(call.rate.batch, cbind(sample = names(cr),
                      call.rate = cr, batch = in.path[i]))
  } # for(i)
  close(x_file)
  close(y_file)
  close(g_file)
  colnames(call.rate.batch) = c("sample", "call.rate", "batch")
  write.table(call.rate.batch, paste(out.path, "call.rate.batch.txt", sep = "/"),
              sep = "\t", row.names = FALSE)

} # extract.raw.data()
