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

  # HDF5 file
  h5tempfile = paste0(tempfile(), ".h5")
  if(!h5createFile(h5tempfile)) {

    stop(paste0("The file \'", h5tempfile, "\' already exists. Please ",
         "move or delete it and try again."))

  } # if(!h5createFile(h5tempfile))

  # Write out the markers.
  h5write(obj = as.matrix(snps), file = h5tempfile, name = "markers")

  # Write out headers for the files.  This will overwrite existing files.
#  x_file = file(paste(out.path, "x.txt", sep = "/"), open = "w")
#  y_file = file(paste(out.path, "y.txt", sep = "/"), open = "w")
#  g_file = file(paste(out.path, "geno.txt", sep = "/"), open = "w")
#  writeLines(text = snps[-nrow(snps),1], con = x_file, sep = "\t")
#  writeLines(text = snps[nrow(snps),1],  con = x_file, sep = "\n")
#  writeLines(text = snps[-nrow(snps),1], con = y_file, sep = "\t")
#  writeLines(text = snps[nrow(snps),1],  con = y_file, sep = "\n")
#  writeLines(text = snps[-nrow(snps),1], con = g_file, sep = "\t")
#  writeLines(text = snps[nrow(snps),1],  con = g_file, sep = "\n")

  call.rate.batch = NULL

  # This will hold a vector of all unique sample IDs for all batches.
  all.samples = NULL
  # Index of the current sample in the all.samples vector.
  index = 1

  for(i in 1:length(in.path)) {

    # Get the sample IDs from the Sample_Map.txt file.
    samplefile = dir(path = in.path[i], pattern = "Sample_Map.txt", full.names = TRUE)

    # If not found, then quit.
    if(length(samplefile) == 0) {
      stop(paste("No file called 'Sample_Map.txt' was found in directory",
           in.path[i], ".  Please make sure that the Sample_Map file is unzipped and",
           "in the specified directory."))
    } # if(length(samplefile) == 0)

    # Read in the sample IDs.
    samples = read.delim(samplefile, stringsAsFactors = FALSE)$Name
    samples = samples[nchar(samples) > 0]
    if(!missing(prefix)) {
      all.samples = make.unique(make.names(c(all.samples, 
                    paste0(prefix[i], samples))))
    } else {
      all.samples = make.unique(make.names(c(all.samples, samples)))
    } # else

    # Find a file with "FinalReport" in the filename.
    reportfile = dir(path = in.path[i], pattern = "FinalReport", full.names = TRUE)
    reportfile = reportfile[grep("txt", reportfile)]

    # If not found, then quit.
    if(length(reportfile) == 0) {
      stop(paste("No file with 'FinalReport' in the filename was found in directory",
           in.path[i], ".  Please make sure that the FinalReport file is unzipped and",
           "in the specified directory."))
    } # if(length(reportfile) == 0)

    # If there is more than one FinalReport file, then quit.
    if(length(reportfile) > 1) {
      stop(paste("There is more than one file with FinalReport in the filename.",
           "Please place only one data set in each directory."))
    } # if(length(reportfile) > 1)

    # Read in the first sample.  The current format requires us to skip 9 lines
    # because there is no comment delimiter at the top of the file.
    print(paste("Reading", reportfile, "..."))
    rawfile = file(reportfile, open = "r")
    data = readLines(con = rawfile, n = 10)
    hdr = strsplit(data, split = "\t")

    # Print number of SNPs and samples.
    print(paste(hdr[[5]], collapse = " "))
    print(paste(hdr[[7]], collapse = " "))
    cn = hdr[[10]]

    # Get number of SNPs.
    num.snps = as.numeric(hdr[[5]][2])

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

    # We read the files in and write them out to conserve memory. As computers
    # get larger, we may be able to keep everything in memory.
    for(j in 1:length(samples)) {

      # Read in the data for one sample.
      data = readLines(con = rawfile, n = num.snps)
      data = strsplit(data, split = "\t")
      data = matrix(unlist(data), nrow = length(data[[1]]), ncol = length(data))
      dimnames(data) = list(cn, data[1,])

      # Check that we have only one sample in this chunk.
      unique.sample = unique(data[rownames(data) == "Sample ID",])
      if(length(unique.sample) > 1) {

        stop(paste("The number of SNPs per sample dose not match the expected",
             "number in", reportfile))

      } # if(length(unique.sample) > 1)

      # Verify that the sample IDs match.
      m = regexpr(samples[j], unique.sample)
      stopifnot(m[1] == 1)

      print(paste("Sample", j, "of", length(samples), ":", samples[j]))

      # Sort the data to match the SNP order.
      data = data[,snps[,1]]

      # X
#      writeLines(samples.in.data[j], con = x_file, sep = "\t")
      xint = as.numeric(data[rownames(data) == "X",])
#      writeLines(xint[-length(xint)], con = x_file, sep = "\t")
#      writeLines(xint[length(xint)], con = x_file)

      # Y
#      writeLines(samples.in.data[j], con = y_file, sep = "\t")
      yint = as.numeric(data[rownames(data) == "Y",])
#      writeLines(yint[-length(yint)], con = y_file, sep = "\t")
#      writeLines(yint[length(yint)], con = y_file)

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
      geno[geno == "--"] = "N"
      geno[geno %in% hets] = "H"

#      writeLines(samples.in.data[j], con = g_file, sep = "\t")
#      writeLines(geno[-length(geno)], con = g_file, sep = "\t")
#      writeLines(geno[length(geno)], con = g_file)

      h5createGroup(file = h5tempfile, group = all.samples[index])
      xy = matrix(c(xint, yint), nrow = length(xint))
      h5write(obj = xy, file = h5tempfile, name = 
              paste0("/", all.samples[index], "/xy"))
      h5write(obj = geno, file = h5tempfile, name = 
              paste0("/", all.samples[index], "/g"))

      index = index + 1

    } # for(j)

    close(rawfile)
    names(cr) = samples
    call.rate.batch = rbind(call.rate.batch, cbind(sample = names(cr),
                      call.rate = cr, batch = in.path[i]))

  } # for(i)

#  close(x_file)
#  close(y_file)
#  close(g_file)

  call.rate.batch[,2] = make.unique(make.names(call.rate.batch[,2]))
  colnames(call.rate.batch) = c("sample", "call.rate", "batch")
  write.table(call.rate.batch, paste(out.path, "call.rate.batch.txt", sep = "/"),
              sep = "\t", row.names = FALSE)

  # Read the HDF5 file back in and create X, Y and Geno matrices.
  # Save these as *.Rdata files.
  info = h5ls(h5tempfile)
  info = info[info$otype == "H5I_GROUP",]

  # Get the unique chromosomes.
  chr = unique(snps[,2])
  chr.factor = factor(snps[,2], levels = (chr), exclude = NULL)
  num.markers.on.chr = sapply(split(snps, chr.factor), nrow)

  # Create X, Y and G matrices for each chromsome.
  # TBD: Wish that I didn't have to read this all into memory.
  x = vector("list", length(chr))
  y = vector("list", length(x))
  g = vector("list", length(x))
  names(x) = chr
  names(y) = names(x)
  names(g) = names(x)

  for(i in 1:length(chr)) {

    x[[i]] = matrix(0, nrow = nrow(info), ncol = num.markers.on.chr[i],
             dimnames = list(info$name, snps[which(snps[,2] == chr[i]),1]))
    y[[i]] = matrix(0, nrow = nrow(x[[i]]), ncol = ncol(x[[i]]), 
             dimnames = dimnames(x[[i]]))
    g[[i]] = matrix("", nrow = nrow(x[[i]]), ncol = ncol(x[[i]]),
             dimnames = dimnames(x[[i]]))

  } # for(i)

  # Read in each sample and place it's data in the correct spot.
  for(i in 1:nrow(info)) {

    xy = data.frame(h5read(file = h5tempfile, name = paste0(info$name[i], "/xy")))
    g2 = h5read(file = h5tempfile, name = paste0(info$name[i], "/g"))

    xy = split(xy, chr.factor)
    g2  = split(g2, chr.factor)

    for(j in 1:length(chr)) {

      x[[j]][i,] = xy[[j]][,1]
      y[[j]][i,] = xy[[j]][,2]
      g[[j]][i,] = g2[[j]]

    } # for(j)

  } # for(i)

  # Create the datasets in the new file.
  h5filename = paste0(out.path, "/x_y_geno.h5")
  if(!h5createFile(h5filename)) {

    stop(paste0("The file \'", h5filename, "\' already exists. Please ",
         "move or delete it and try again."))

  } # if(!h5createFile(h5filename))

  h5file = H5Fopen(name = h5filename)

  # Write out the markers.
  h5write(obj = as.matrix(snps), file = h5filename, name = "markers")

  # Write out the sample IDs.
  h5write(obj = info$name, file = h5filename, name = "samples")

  # Write out the call rate and batch information.
  h5write(obj = as.matrix(call.rate.batch), file = h5filename, name = "batch")

  # Create the data sets, with chunking by samples if there are more than
  # 100 samples.
  h5createGroup(file = h5filename, group = "X")
  h5createGroup(file = h5filename, group = "Y")
  h5createGroup(file = h5filename, group = "G")

  # Set the chunk size to 100 samples if num.samples > 100)
  chunk = nrow(info)
  if(nrow(info) > 100) {
   chunk = 100
  } # if(nrow(info) > 100)

  # Create a dataset for each chromosome in X, Y and G.
  for(i in chr[!is.na(chr)]) {

    print(paste("Writing Chr", i))

    # X
    h5createDataset(file = h5filename, dataset = paste0("/X/", i),
                    dims = c(nrow(info), num.markers.on.chr[i]),
                    chunk = c(chunk, num.markers.on.chr[i]), showWarnings = FALSE)
    h5grp = H5Gopen(h5loc = h5file, name = "/X")
    h5writeDataset(obj = x[[i]], h5loc = h5grp, name = i)

    # Y
    h5createDataset(file = h5filename, dataset = paste0("/Y/", i),
                    dims = c(nrow(info), num.markers.on.chr[i]),
                    chunk = c(chunk, num.markers.on.chr[i]), showWarnings = FALSE)
    h5grp = H5Gopen(h5loc = h5file, name = "/Y")
    h5writeDataset(obj = y[[i]], h5loc = h5grp, name = i)

    # G
    h5createDataset(file = h5filename, dataset = paste0("/G/", i),
                    dims = c(nrow(info), num.markers.on.chr[i]),
                    storage.mode = "character", size = 1,
                    chunk = c(chunk, num.markers.on.chr[i]), showWarnings = FALSE)
    h5grp = H5Gopen(h5loc = h5file, name = "/G")
    h5writeDataset(obj = g[[i]], h5loc = h5grp, name = i)

  } # for(i)

  H5close()

} # extract.raw.data()
