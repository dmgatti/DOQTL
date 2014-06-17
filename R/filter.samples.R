################################################################################
# Filter the samples by call rate.
# The function attempts to figure out if the user is working with intensity 
# or genotype data by looking for the x.txt and y.txt files.
# Arguments: path: character vector with full path to data files.
#            thr: numeric, call rate threshold.
# Returns: data.frame of samples removed.
# Output: genotypes.filt.txt and call.rate.batch.filt.txt with samples removed
#         with call rates below the threshold. Possibly x.filt.txt and
#         y.filt.txt if x.txt and y.txt were found in the path.
################################################################################
filter.samples <- function(path = ".", thr = 0.9) {
  print("Reading files..")
  crb = read.delim(paste(path, "call.rate.batch.txt", sep = "/"), sep = "\t")
  xfile = paste(path, "x.txt", sep = "/")
  yfile = paste(path, "y.txt", sep = "/")
  gfile = paste(path, "geno.txt", sep = "/")
  xin = file(xfile, open = "r") 
  yin = file(yfile, open = "r") 
  gin = file(gfile, open = "r") 
  xout = file(paste(path, "x.filt.txt", sep = "/"), open = "w") 
  yout = file(paste(path, "y.filt.txt", sep = "/"), open = "w") 
  gout = file(paste(path, "geno.filt.txt", sep = "/"), open = "w") 
  # Copy the header line to the output files.
  x = readLines(xin, n = 1)
  y = readLines(yin, n = 1)
  g = readLines(gin, n = 1)
  writeLines(x, con = xout)  
  writeLines(y, con = yout)
  writeLines(g, con = gout)
  print("Writing files..")
  for(i in 1:nrow(crb)) {
    print(paste(i, "of", nrow(crb)))
    x = readLines(xin, n = 1)
    y = readLines(yin, n = 1)
    g = readLines(gin, n = 1)
    tmp = strsplit(x, split = "\t")
    if(tmp[[1]][1] != crb[i,1]) {
      close(xin); close(yin); close(gin)
      close(xout); close(yout); close(gout)
      stop(paste("Sample IDs are mismatched between the call rate file",
           "and the x, y and g files."))
    } # if(tmp[[1]][1] != crb[i,1])
    
    if(crb[i,2] >= thr) {
      writeLines(x, con = xout)  
      writeLines(y, con = yout)
      writeLines(g, con = gout)
    } # if(crb[i,2] >= thr)
  } # for(i)
  close(xin)
  close(yin)
  close(gin)
  close(xout)
  close(yout)
  close(gout)
  retval = crb[crb[,2] < thr,]
  crb = crb[crb[,2] >= thr,]
  # Write out the call rate/batch file.
  write.table(crb, paste(path, "call.rate.batch.filt.txt", sep = "/"), 
              sep = "\t")
  return(retval)
} # filter.samples()
