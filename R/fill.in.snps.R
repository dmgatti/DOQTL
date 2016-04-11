################################################################################
# Go through the SNPs and look for stretches where the differences in cM values
# from one SNP to the next == 0. Interpolate from the first SNP to the last SNP
# in such a stretch.
# Arguments:
# snps: data.frame with 4 columns. Column 1: SNP IDs, column 2: chr,
#          column 3: Mb position, column 4: cM position.
fill.in.snps = function(snps) {

  # Order the SNPs.
  num.auto = max(as.numeric(snps[,2]), na.rm = TRUE)
  snps[snps[,2] == "X",2] = num.auto + 1
  snps[snps[,2] == "Y",2] = num.auto + 2
  snps[snps[,2] == "M",2] = num.auto + 3
  snps = snps[order(snps[,3]),]
  snps[,2] = as.numeric(snps[,2])
  snps = snps[order(snps[,2]),]
  # Find the SNPs where the cM do not differ from one SNP to the next.  
  zero = c(diff(snps[,4]) == 0, 0)
  zero = c(0, diff(zero))
  start = which(zero ==  1)
  end   = which(zero == -1)

  if(snps[1,4] == snps[2,4]) {
    start = c(1, start)
  } # if(snps[1,4] == snps[2,4])
  
  if(length(start) > 0 & length(end) > 0) {
    # Find the Chr boundaries.
    chr.bounds = table(snps[,2])
    old.warn = options("warn")$warn
    options(warn = FALSE)
    chr.bounds = chr.bounds[order(as.numeric(names(chr.bounds)))]
    options(warn = old.warn)
    chr.bounds = cumsum(chr.bounds)
    chr.bounds = c(1, chr.bounds, chr.bounds + 1)
    # Handle the Chr ends first.
    st = which(start %in% chr.bounds)
    if(length(st) > 0) {
      for(i in st) {
        step = (snps[end[i] + 1,4] - snps[start[i],4]) / 
               (end[i] - start[i] + 1)
        snps[(start[i] + 1):(end[i]), 4] = snps[start[i],4] + 
                                           1:(end[i] - start[i]) * step
      } # for(i)
    } # if(length(st) > 0)
    ed = which(end %in% chr.bounds)
    if(length(ed) > 0) {
      for(i in ed) {
        step = (snps[end[i],4] - snps[start[i] - 1,4]) / 
               (end[i] - start[i] + 1)
        snps[(start[i]):(end[i] - 1), 4] = snps[start[i] - 1,4] + 
                                           1:(end[i] - start[i]) * step
      } # for(i)
    } # if(length(st) > 0)
    if(length(ed) > 0) {
      start = start[-ed]
      end   = end[-ed]
    } # if(length(ed) > 0)
    if(length(st) > 0) {
      start = start[-st]
      end   = end[-st]
    } # if(length(st) > 0)
    if(length(start) > 0) {
      for(i in 1:length(start)) {
        step = (snps[end[i] + 1,4] - snps[start[i] - 1,4]) / 
               (end[i] - start[i] + 2)
        snps[start[i]:end[i], 4] = snps[start[i] - 1,4] + 
                                   1:(end[i] - start[i] + 1) * step
      } # for(i)
    } # if(length(start) > 0)
  } # if(length(start) > 0 & length(end) > 0)
  snps[snps[,2] == num.auto + 1,2] = "X"
  snps[snps[,2] == num.auto + 2,2] = "Y"
  snps[snps[,2] == num.auto + 3,2] = "M"
  return(snps)
} # fill.in.snps()
