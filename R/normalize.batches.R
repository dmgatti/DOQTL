# Normalize batches by centering on median overall intensity.
# Arguments: x: matrix of X intensities, samples in rows, SNPs in columns.
#            y: matrix of Y intensities, samples in rows, SNPs in columns.
#            sex: vector of "F" or "M" for each sample.  Required for X and 
#                 Y Chr normalization.
#            batch: vector, factor of batch identifiers.
#            snps: data.frame with three columns: SNP ID, chr, location.
normalize.batches = function(x, y, sex, batch, snps) {

  batch = factor(batch)
  sample.size.cutoff = 30

  # Scale all of the data to be between 0 and 1.
  x = t(x)
  x = x - apply(x, 1, min, na.rm = T)
  x = x / apply(x, 1, max, na.rm = T)
  x = t(x)
  y = t(y)
  y = y - apply(y, 1, min, na.rm = T)
  y = y / apply(y, 1, max, na.rm = T)
  y = t(y)

  # First adjust the autosomes.
  old.warn = options("warn")$warn
  options(warn = -1)
  auto.snps = which(!is.na(as.numeric(snps[,2])))
  options(warn = old.warn)

  # Compute the overall intensity as X + Y intensity.
  intensity = x + y
  overall.med = apply(intensity, 2, median, na.rm = T)

  for(i in levels(batch)) {

    rng = which(batch == i)
    if(length(rng) > sample.size.cutoff) {
      batch.factor = overall.med[auto.snps] /
                     apply(intensity[rng, auto.snps, drop=F], 2, median,
                     na.rm = T)
      batch.factor = matrix(batch.factor, length(rng), length(auto.snps),
                            byrow = T)
      x[rng, auto.snps] = x[rng, auto.snps, drop=F] * batch.factor
      y[rng, auto.snps] = y[rng, auto.snps, drop=F] * batch.factor
    } # if(length(rng) > sample.size.cutoff)
  } # for(i)

  # Next, adjust the males and females separately on the X & Y chromosomes.
  x.snps = which(snps[,2] == "X")
  y.snps = which(snps[,2] == "Y")

  female.samples = which(sex == "F")
  male.samples   = which(sex == "M")

  if(length(female.samples) > 0) {
    f.med = apply(intensity[female.samples,, drop=F], 2, median, na.rm = T)

    for(i in levels(batch)) {
      # Select the female samples in the current batch.
      rng = intersect(which(batch == i), female.samples)
      if(length(rng) > sample.size.cutoff) {
        # Create a factor by which to multiply the X & Y intensities.
        batch.factor = f.med[x.snps] /
                       apply(intensity[rng, x.snps, drop=F], 2, median,
                       na.rm = T)
        batch.factor = matrix(batch.factor, length(rng), length(x.snps),
                       byrow = T)
        x[rng, x.snps] = x[rng, x.snps, drop=F] * batch.factor
        y[rng, x.snps] = y[rng, x.snps, drop=F] * batch.factor

        batch.factor = f.med[y.snps] /
                       apply(intensity[rng, y.snps, drop=F], 2, median,
                       na.rm = T)
        batch.factor = matrix(batch.factor, length(rng), length(y.snps),
                       byrow = T)
        x[rng, y.snps] = x[rng, y.snps, drop=F] * batch.factor
        y[rng, y.snps] = y[rng, y.snps, drop=F] * batch.factor
      } # if(length(rng) > sample.size.cutoff)
    } # for(i)
  } # if(length(female.samples) > 0)

  if(length(male.samples) > 0) {
    m.med = apply(intensity[male.samples, , drop=F], 2, median, na.rm = T)

    for(i in levels(batch)) {
      # Select the male samples in the current batch.
      rng = intersect(which(batch == i), male.samples)
      if(length(rng) > sample.size.cutoff) {
        # Create a factor by which to multiply the X & Y intensities.
        batch.factor = m.med[x.snps] /
                       apply(intensity[rng, x.snps, drop=F], 2, median,
                       na.rm = T)
        batch.factor = matrix(batch.factor, length(rng), length(x.snps),
                              byrow = T)
        x[rng, x.snps] = x[rng, x.snps, drop=F] * batch.factor
        y[rng, x.snps] = y[rng, x.snps, drop=F] * batch.factor

        batch.factor = m.med[y.snps] /
                       apply(intensity[rng, y.snps, drop=F], 2, median,
                       na.rm = T)
        batch.factor = matrix(batch.factor, length(rng), length(y.snps),
                              byrow = T)
        x[rng, y.snps] = x[rng, y.snps, drop=F] * batch.factor
        y[rng, y.snps] = y[rng, y.snps, drop=F] * batch.factor
      } # if(length(rng) > sample.size.cutoff)
    } # for(i)
  } # if(length(male.samples) > 0)

  return(list(x = x, y = y))
} # normalize.batches()


# Helper function
batch.normalize = function(path = ".", snps) {

  # Figure out if we have filtered samples or not.
  # X data
  x.file = dir(path, pattern = "^x.txt", full.names = T)
  x.filt.file = dir(path, pattern = "^x.filt.txt", full.names = T)
  if(length(x.filt.file) > 0) {
    x.file = x.filt.file
  } # if(length(x.filt.file) > 0)

  x = read.delim(x.file, row.names = NULL)
  rn = x[,1]
  x = as.matrix(x[,-1])
  rownames(x) = rn

  # Y data
  y.file = dir(path, pattern = "y.txt", full.names = T)
  y.filt.file = dir(path, pattern = "y.filt.txt", full.names = T)
  if(length(y.filt.file) > 0) {
    y.file = y.filt.file
  } # if(length(y.filt.file) > 0)

  y = read.delim(y.file, row.names = NULL)
  rn = y[,1]
  y = as.matrix(y[,-1])
  rownames(y) = rn

  # Batch data
  crb.file = dir(path, pattern = "call.rate.batch.txt", full.names = T)
  crb.filt.file = dir(path, pattern = "call.rate.batch.filt.txt", full.names = T)
  if(length(crb.filt.file) > 0) {
    crb.file = crb.filt.file
  } # if(length(crb.filt.file) > 0)

  crb = read.delim(crb.file)
  crb$batch = factor(crb$batch)

  if(nrow(x) != nrow(y)) {
    stop(paste("The number of rows in X is not equal to the number of rows in Y."))
  } # if(nrow(x) != nrow(y))

  if(ncol(x) != ncol(y)) {
    stop(paste("The number of columns in X is not equal to the number of columns in Y."))
  } # if(ncol(x) != ncol(y))

  if(nrow(x) != nrow(crb)) {
    stop(paste("The number of rows in X is not equal to the number of rows in",
         "the call rates."))
  } # if(nrow(x) != nrow(y))

  # Estimate the sex of each sample.
  sex = predict.sex(x, y, snps, plot = F)

  # Compute the overall intensity as X + Y intensity.
  tmp = normalize.batches(x, y, sex, crb$batch, snps) 

  write.table(tmp$x, sub("txt$", "batch.norm.txt", x.file), sep = "\t")
  write.table(tmp$y, sub("txt$", "batch.norm.txt", y.file), sep = "\t")
} # batch.normalize()


