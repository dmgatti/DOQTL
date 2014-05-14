# Normalize batches by centering on 90th %-tile intensity.
# Arguments: x: matrix of X intensities, samples in rows, SNPs in columns.
#            y: matrix of Y intensities, samples in rows, SNPs in columns.
#            sex: vector of "F" or "M" for each sample.  Required for X and 
#                 Y Chr normalization.
#            batch: vector, factor of batch identifiers.
#            snps: data.frame with three columns: SNP ID, chr, location.
normalize.batches = function(x, y, sex, batch, snps) {
  batch = factor(batch)
  sample.size.cutoff = 20
  # First adjust the autosomes.
  old.warn = options("warn")$warn
  options(warn = -1)
  auto.snps = which(!is.na(as.numeric(snps[,2])))
  options(warn = old.warn)
  # Compute the overall X and Y 90th %-tiles for each autosomal SNP.
  x.q9 = apply(x[,auto.snps], 2, quantile, probs = 0.9, na.rm = TRUE)
  y.q9 = apply(y[,auto.snps], 2, quantile, probs = 0.9, na.rm = TRUE)
  for(i in levels(batch)) {
    rng = which(batch == i)
    if(length(rng) > sample.size.cutoff) {
       x[rng, auto.snps] = x[rng, auto.snps, drop = FALSE] / 
                 matrix(apply(x[rng, auto.snps, drop = FALSE], 2, quantile, 
                 0.9, na.rm = TRUE), length(rng), length(auto.snps), byrow = TRUE) *
                 x.q9
       y[rng, auto.snps] = y[rng, auto.snps] / 
                 matrix(apply(y[rng, auto.snps, drop = FALSE], 2, quantile,
                 0.9, na.rm = TRUE), length(rng), length(auto.snps), byrow = TRUE) * 
                 y.q9
    } # if(length(rng) > sample.size.cutoff)
  } # for(i)
  # Next, adjust the males and females separately on the X & Y chromosomes.
  x.snps = which(snps[,2] == "X")
  y.snps = which(snps[,2] == "Y")
  female.samples = which(sex == "F")
  male.samples   = which(sex == "M")
  intensity = x + y
  
  if(length(female.samples) > 0) {
    f.med = apply(intensity[female.samples,, drop=FALSE], 2, median, na.rm = TRUE)
    for(i in levels(batch)) {
      # Select the female samples in the current batch.
      rng = intersect(which(batch == i), female.samples)
      if(length(rng) > sample.size.cutoff) {
        # Create a factor by which to multiply the X & Y intensities.
        batch.factor = f.med[x.snps] /
                       apply(intensity[rng, x.snps, drop=FALSE], 2, median,
                       na.rm = TRUE)
        batch.factor = matrix(batch.factor, length(rng), length(x.snps),
                       byrow = TRUE)
        x[rng, x.snps] = x[rng, x.snps, drop=FALSE] * batch.factor
        y[rng, x.snps] = y[rng, x.snps, drop=FALSE] * batch.factor
        batch.factor = f.med[y.snps] /
                       apply(intensity[rng, y.snps, drop=FALSE], 2, median,
                       na.rm = TRUE)
        batch.factor = matrix(batch.factor, length(rng), length(y.snps),
                       byrow = TRUE)
        x[rng, y.snps] = x[rng, y.snps, drop=FALSE] * batch.factor
        y[rng, y.snps] = y[rng, y.snps, drop=FALSE] * batch.factor
      } # if(length(rng) > sample.size.cutoff)
    } # for(i)
  } # if(length(female.samples) > 0)
  if(length(male.samples) > 0) {
    m.med = apply(intensity[male.samples, , drop=FALSE], 2, median, na.rm = TRUE)
    for(i in levels(batch)) {
      # Select the male samples in the current batch.
      rng = intersect(which(batch == i), male.samples)
      if(length(rng) > sample.size.cutoff) {
        # Create a factor by which to multiply the X & Y intensities.
        batch.factor = m.med[x.snps] /
                       apply(intensity[rng, x.snps, drop=FALSE], 2, median,
                       na.rm = TRUE)
        batch.factor = matrix(batch.factor, length(rng), length(x.snps),
                              byrow = TRUE)
        x[rng, x.snps] = x[rng, x.snps, drop=FALSE] * batch.factor
        y[rng, x.snps] = y[rng, x.snps, drop=FALSE] * batch.factor
        batch.factor = m.med[y.snps] /
                       apply(intensity[rng, y.snps, drop=FALSE], 2, median,
                       na.rm = TRUE)
        batch.factor = matrix(batch.factor, length(rng), length(y.snps),
                              byrow = TRUE)
        x[rng, y.snps] = x[rng, y.snps, drop=FALSE] * batch.factor
        y[rng, y.snps] = y[rng, y.snps, drop=FALSE] * batch.factor
      } # if(length(rng) > sample.size.cutoff)
    } # for(i)
  } # if(length(male.samples) > 0)
  return(list(x = x, y = y))
} # normalize.batches()
# Helper function that reads in the files.
batch.normalize = function(path = ".", snps) {
  message("Reading in data ...")
  # Figure out if we have filtered samples or not.
  # X data
  x.file = dir(path, pattern = "^x.txt", full.names = TRUE)
  x.filt.file = dir(path, pattern = "^x.filt.txt", full.names = TRUE)
  if(length(x.filt.file) > 0) {
    x.file = x.filt.file
  } # if(length(x.filt.file) > 0)
  x = read.delim(x.file, row.names = NULL)
  rn = x[,1]
  x = as.matrix(x[,-1])
  rownames(x) = rn
  # Y data
  y.file = dir(path, pattern = "y.txt", full.names = TRUE)
  y.filt.file = dir(path, pattern = "y.filt.txt", full.names = TRUE)
  if(length(y.filt.file) > 0) {
    y.file = y.filt.file
  } # if(length(y.filt.file) > 0)
  y = read.delim(y.file, row.names = NULL)
  rn = y[,1]
  y = as.matrix(y[,-1])
  rownames(y) = rn
  # Batch data
  crb.file = dir(path, pattern = "call.rate.batch.txt", full.names = TRUE)
  crb.filt.file = dir(path, pattern = "call.rate.batch.filt.txt", full.names = TRUE)
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
  sex = sex.predict(x, y, snps, plot = FALSE)
  # Order the batches by sample size, largest to smallest.
  batch.order = sort(table(crb$batch), decreasing = TRUE)
  # Remove batches with less than 50 samples because we can't reliably
  # normalize small batches. We will leave small batches unchanged.
  batch.order = batch.order[batch.order >= 45]
  message("Normalizing batches ...")
  # Normalize to the largest batch.
  for(i in 2:length(batch.order)) {
    b1 = which(crb$batch %in% names(batch.order)[1])
    b2 = which(crb$batch == names(batch.order)[i])
    newxy = quantilenorm(x1 = x[b1,], y1 = y[b1,], x2 = x[b2,], y2 = y[b2,])
    x[b2,] = newxy[[1]]
    y[b2,] = newxy[[2]]
  } # for(i)
  write.table(x, sub("txt$", "batch.norm.txt", x.file), sep = "\t")
  write.table(y, sub("txt$", "batch.norm.txt", y.file), sep = "\t")
} # batch.normalize()
##########
# This function uses a quantile normalization separately for the X and Y
# intensities at each SNP.
# We assume that data set 1 (x1 & y1) contain more samples than data set 2
# (x2 & y2) and adjust data set 2 to data set 1.
# Arguments: x1: numeric matrix containing X intensities for batch 1
#                containing samples in rows and markers in columns.
#            y1: numeric matrix containing Y intensities for batch 1 
#                containing samples in rows and markers in columns.
#            x2: numeric matrix containing X intensities for batch 2
#                containing samples in rows and markers in columns.
#            y2: numeric matrix containing Y intensities for batch 2
#                containing samples in rows and markers in columns.
quantilenorm = function(x1, y1, x2, y2) {
  # Error checking.
  if(ncol(x1) != ncol(y1)) {
    stop(paste("The number of columns in x1 does not equal the",
         "number columns in y1."))
  } # if(ncol(x1) != ncol(y1))
  if(ncol(x1) != ncol(x2)) {
      stop(paste("The number of columns in x1 does not equal the",
         "number columns in x2."))
  } # if(ncol(x1) != ncol(x2))
  if(ncol(y1) != ncol(y2)) {
      stop(paste("The number olf columns in y1 does not equal the",
         "number columns in y2."))
  } # if(ncol(y1) != ncol(y2))
  # Convert x1 & y1 matrices to data.frames.
  if(!is.data.frame(x1)) {
    x1 = as.data.frame(x1)
    y1 = as.data.frame(y1)
  } # if(!is.matrix(x1))
  num.snps = ncol(x1)
  # We need to remove NAs here because we need real values at each quantile.
  x1 = lapply(x1, sort, na.last = NA)
  y1 = lapply(y1, sort, na.last = NA)
  # At this point, x1 and y1 are lists.
  # Normalize at each SNP.
  for(s in 1:num.snps) {
    if(s %% 1000 == 0) print(paste(s/num.snps  *100, "%"))
  
    # Keep only samples between the 0.01 and 0.99 quantiles to reduce the
    # effect of a handful of outliers.
    keep = ceiling(0.01 * length(x1[[s]])):floor(0.99 * length(x1[[s]]))
    # Prepare X1 & Y1 by aligning quantiles and intensity values.
    sortval = x1[[s]][keep]
    quantval = sortval - sortval[1]
    qx1 = cbind(quantval / quantval[length(quantval)], sortval)
    qx1 = qx1[!duplicated(qx1[,1]),,drop=FALSE]
    # This occurs when all values of x are the same (often 0).
    if(nrow(qx1) == 1) {
      qx1 = matrix(c(0, 1, rep(sortval[1], 2)), 2, 2)
    } # if(nrow(qx1) == 1)
    sortval = y1[[s]][keep]
    quantval = sortval - sortval[1]
    qy1 = cbind(quantval / quantval[length(quantval)], sortval)
    qy1 = qy1[!duplicated(qy1[,1]),,drop=FALSE]
    # This occurs when all values of y are the same (often 0).
    if(nrow(qy1) == 1) {
      qy1 = matrix(c(0, 1, rep(sortval[1], 2)), 2, 2)
    } # if(nrow(qx1) == 1)
    
    # Prepare X2 & Y2 by aligning quantiles, intensity values, new
    # intensity values and the original sort order.
    ord = order(x2[,s], na.last = FALSE)
    sortval = x2[ord,s]
    quantval = sortval - min(sortval, na.rm = TRUE)
    qx2 = cbind(quantval / max(quantval, na.rm = TRUE), sortval, 
          rep(0, length(sortval)), ord)
    ord = order(y2[,s], na.last = FALSE)
    sortval = y2[ord,s]
    quantval = sortval - min(sortval, na.rm = TRUE)
    qy2 = cbind(quantval / max(quantval, na.rm = TRUE), sortval, 
          rep(0, length(sortval)), ord)
    # Adjust X by finding the quantile in x1 where each x2 point occurs.
    # lt is the index of the quantile less than the x2 point.
    # gt is the index of the quantile greater than the x2 point.
    diff = outer(qx1[,1], qx2[,1], "-")
    lt = colSums(diff <= 0)
    lt[lt < 1] = 1
    lt[lt == nrow(qx1)] = nrow(qx1) - 1
    gt = lt + 1
    qx2[,3] = qx1[lt,2] + (qx1[gt,2] - qx1[lt,2]) * 
              (qx2[,1] - qx1[lt,1]) / (qx1[gt,1] - qx1[lt,1])
    x2[qx2[,4],s] = qx2[,3]
    # Adjust Y by finding the quantile in y1 where each y2 point occurs.
    diff = outer(qy1[,1], qy2[,1], "-")
    lt = colSums(diff < 0)
    lt[lt < 1] = 1
    lt[lt == nrow(qy1)] = nrow(qy1) - 1
    gt = lt + 1
    qy2[,3] = qy1[lt,2] + (qy1[gt,2] - qy1[lt,2]) * 
              (qy1[gt,1] - qy2[,1]) / (qy1[gt,1] - qy1[lt,1])
    y2[qy2[,4],s] = qy2[,3]
  } # for(s)
  return(list(x2, y2))
} # quantilenorm()
