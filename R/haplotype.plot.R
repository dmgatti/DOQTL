################################################################################
# Given the genotype probabilities and a set of alleles, plot a region where
# all samples have the requested alleles. This might be used to narrow a QTL
# interval where a dominant effect is suspected.
# Daniel Gatti
# Dan.Gatti@jax.org
# Aug. 4, 2013
################################################################################
# Arguments: probs: numeric 3D array containing the founder probabilities with
#                   dimensions samples x founders x SNPs. All dimensions must
#                   have dimnames.
#            snps: data.frame containing SNP information. SNP ID, Chr, Mb in 
#                  columns 1:3.  SNP IDs much match those in dimnames(probs)[[3]].
#            chr: character containing a single chromosome.
#            start: numeric value containing the starting postiton for plotting
#                   in Mb.
#            end: numeric value containing the ending postiton for plotting
#                 in Mb.
#            strain: character contianing a single letter code that is in 
#                    dimnames(probe)[[2]] that is the strain haplotype to plot.
#            colors: either "DO", in which case the DO founder colors will be
#                    used, or a data.frame with founder letters, founder names,
#                    founder colors in columns 1:3.
#            ...: Additional arguments to be passed into plot.
haplotype.plot = function(probs, snps, chr, start, end, strain, colors = "DO",
                          ...) {
  
  if(length(chr) > 1) {
    warning("Using only the first element of chr.")
    chr = chr[1]
  } # if(length(chr) > 1)
  
  if(length(strain) > 1) {
    warning(paste("Using only the first element of strain. For more strains,",
            "please call the function once per strain."))
  } # if(length(strain) > 1)
  
  # Select the desired interval
  snps = snps[snps[,2] == chr,]
  snps = snps[snps[,3] >= start & snps[,3] <= end,]
  
  # Subset the probs and SNPS.
  probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
  
  # Keep only the deisred strains.
  probs = probs[,dimnames(probs)[[2]] %in% strain,]
  
  if(colors == "DO") {
    data(founder.cols)
    colors = founder.cols
  } else {
    if(!is.data.frame(colors)) {
      stop(paste("If not using the DO colors, colors must be a data frame",
                 "with three columns containing the founder letters, names and",
                 "colors in columns 1, 2 & 3, respectively."))
    } # if(!is.data.frame(colors))
    
    if(ncol(colors) != 3) {
      stop(paste("Colors does not have three columns. Colors must have three",
                 "columns containing the founder letters, names and colors in",
                 "columns 1, 2 & 3, respectively."))
    } # if(ncol(colors) != 3)
  } # else
  
  # Sort the samples by their haplotypes.
  probs.aggr = apply(probs, 1, sum)
  probs = probs[order(probs.aggr, decreasing = F),]
  

  par(las = 1)
  breaks = c(-0.25, 0.25, 0.75, 1.25)
  col = colorRampPalette(c(colors[colors[,1] == strain,3], rgb(0,0,0)))(length(breaks) - 1)
  image(snps[,3], 1:nrow(probs), t(probs), breaks = breaks, col = col, ylab = "Samples", xlab = 
       paste("Chr", chr, "(Mb)"), main = founder.cols[founder.cols[,1] == strain ,2], ...)
  abline(h = 0:nrow(probs) + 0.5, col = "grey50")
  abline(v = round(snps[,3]), col = "grey50")
} # haplotype.plot()