################################################################################
# Write out the smoothed probabilities for each sample and the founder state
# means and variances. This will write out a file for all samples on the 
# autosomes.  On the X Chr, you must specify which sex is being written out.
# 
write.results = function(prsmth, theta.rho.means, theta.rho.covars, b, output.dir,
                         chr, all.chr, sex) {
  if(!missing(prsmth)) {
    dimnames(prsmth)[[2]] = make.names(dimnames(prsmth)[[2]])
    for(i in 1:dim(prsmth)[2]) {
       write.table(t(prsmth[,i,]), paste(output.dir, dimnames(prsmth)[[2]][i],
                   ".genotype.probs.txt", sep = ""), append = chr != all.chr[1],
                   sep = "\t", col.names = chr == all.chr[1]) 
    } # for(i)
  } # if(!missing(prsmth))
  # Write out the genotype mean and variance estimates.
  if(!missing(theta.rho.means)) {  
    if(chr != "X") {
      # Write out the genotype mean and variance estimates.
      save(theta.rho.means, file = paste(output.dir, "chr", chr,
           ".final.means.Rdata", sep = ""))
      save(theta.rho.covars, file = paste(output.dir, "chr", chr,
           ".final.covars.Rdata", sep = ""))
    } else {
    
      if(is.null(sex)) {
        stop("write.results: sex cannot be null in write.results if chr = X.")
      } # if(is.null(sex))
    
      if(sex == "M") {
        sex = "male"
      } else {
        sex = "female"
      } # else
    
      save(theta.rho.means, file = paste(output.dir, "chr", chr, ".", sex,
           ".final.means.Rdata", sep = ""))
      save(theta.rho.covars, file = paste(output.dir, "chr", chr, ".", sex,
           ".final.covars.Rdata", sep = ""))
    } # else
  } else if(!missing(b)) {
    if(chr != "X") {
      # Write out the genotype mean and variance estimates.
      save(b, file = paste(output.dir, "chr", chr, ".emission.probs.Rdata", 
	       sep = ""))
    } else {
    
      if(is.null(sex)) {
        stop("write.results: sex cannot be null in write.results if chr = X.")
      } # if(is.null(sex))
    
      if(sex == "M") {
        sex = "male"
      } else {
        sex = "female"
      } # else
    
      save(b, file = paste(output.dir, "chr", chr, ".", sex, 
	       ".emission.probs.Rdata", sep = ""))
    } # else    
  } # else if(!missing(b))
  
} # write.results()
