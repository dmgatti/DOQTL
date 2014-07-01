################################################################################
# Read Sanger VCF files.
# Daniel Gatti
# Dan.Gati@jax.org
# Oct. 24, 2013
# NOTE: Function under development. Currently only handles SNP files.
################################################################################
# Arguments: vcf.file: character string with full path to Sanger VCF file.
#            chr: character or numeric vector indicating chromosomes.
#            start: numeric vector containing start positions (in MB or bp).
#            end: numeric vector containing end positions (in MB or bp).
#            strains: chracter vector for strain names to retreive.
#            return.val: character indicating the kind of values to return. 
#                    If "alleles", in which case the allele call letters will be
#                    returned. If "numbers", the numeric values for the alleles
#                    as found in the Sanger files.
#            return.qual: boolean that is TRUE if the quality scores should be
#                         returned.
#            csq: boolean that is TRUE if variant consequence should be returned.
read.vcf = function(vcf.file, chr = 1, start = 4, end = 4.5, strains,
           return.val = c("allele", "number"), return.qual = TRUE, csq = FALSE) {
		   
  return.val = match.arg(return.val)
  
  if(missing(vcf.file)) {
    stop(paste("vcf.file cannot be missing. Please enter the full path to a",
         "Sanger Mouse Genomes Project VCF file."))
  } # if(missing(vcf.file))
  
  if(length(chr) != length(start) | length(chr) != length(end)) {
    stop(paste("The nubmer of locations in chr, start and end must be equal.\n",
         "len(chr) =", length(chr), "len(start) =", length(start), "len(end) =",
         length(end)))
  } # if(length(chr) != length(start) | ...
  
  if(any(start > end)) {
    stop(paste("Start cannot be larger than end. Please enter a start position",
         "that is less than the end position."))
  } # if(start > end)

  # In the mouse, we assume that start values less than 200 are in Mb because
  # the longest chromosome is less than 200 Mb.
  if(any(start <= 200)) {
    start = start * 1e6
  } # if(start <= 200)

  if(any(end <= 200)) {
    end = end * 1e6
  } # if(end <= 200)

  # Query Tabix indexed VCF file.
  tabix = TabixFile(file = vcf.file)
  open(con = tabix)
  hdr = headerTabix(file = tabix)
  gr = GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  data = scanTabix(file = tabix, param = gr)
  close(con = tabix)
  colnames = sub(hdr$comment, "", hdr$header[length(hdr$header)])
  colnames = strsplit(colnames,, split = "\t")[[1]]
  keep = 1:length(colnames)
  
  if(!missing(strains)) {
    if(!all(strains %in% colnames)) {
      stop(paste("All strains not in VCF file. Missing strains:",
           paste(strains[!strains %in% colnames], collapse = ",")))
    } # if(!all(strains %in% colnames))
    keep = c(1:9, which(colnames %in% strains))
  } # if(!missing(strains))
  
  if(length(grep("snp", vcf.file)) > 0) {
    retval = read.vcf.snp(data, keep, colnames, csq, return.val, return.qual)
  } else if(length(grep("indel", vcf.file)) > 0) {
    retval = read.vcf.indel(data, keep, colnames, csq, return.val, return.qual)
  } else if(length(grep("SV", vcf.file)) > 0) {
    retval = read.vcf.sv(data, keep, colnames, csq, return.val, return.qual)
  } else {
    stop(paste("Unknown file type. The VCF file name must contain one of",
         "snp, indel or SV in order to identify the type of file to parse."))
  } # else
  
  # If there was only one query, then return a data.frame rather than a list
  # of data.frames.
  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)
  return(retval)
  
} # read.vcf()
################################################################################
# Retrieve SNPs from a Sanger VCF file.
################################################################################
read.vcf.snp = function(data, keep, colnames, csq, return.val, return.qual) {

  data.not.empty = which(sapply(data, length) > 0)
  
  for(i in data.not.empty) {
  
    data[[i]] = strsplit(data[[i]], split = "\t")
    data[[i]] = matrix(unlist(data[[i]]), nrow = length(data[[i]]),
                byrow = TRUE, dimnames = list(NULL, colnames))
    data[[i]] = data[[i]][,keep,drop = FALSE]
    
    # Split up the genotype calls and quality scores.
    d2 = apply(data[[i]][,10:ncol(data[[i]]),drop=FALSE], 2, strsplit, split = ":")
#    d2 = lapply(d2, function(z) { lapply(z, function(a) { a[c(1,6)] }) })
    d2 = lapply(d2, lapply, "[", c(1,6))
    geno = matrix(unlist(lapply(d2, sapply, "[", 1)),
           ncol = length(d2), dimnames = list(NULL, names(d2)))
    qual = matrix(unlist(lapply(d2, sapply, "[", 2)),
           ncol = length(d2), dimnames = list(NULL, names(d2)))
    rm(d2)
    gc()
	
    if(return.val == "number") {
	
      geno = sub("/", "", geno)
      geno = sub("\\.", NA, geno)
      geno = matrix(as.numeric(geno), nrow(geno), ncol(geno), dimnames = 
             dimnames(geno))
			 
    } else if(return.val == "allele") {
	
      # Create a matrix that maps numeric allele calls to nucleotide allele
      # calls.
      alt = strsplit(data[[i]][,colnames(data[[i]]) == "ALT"], split = ",")
      ref.col = which(colnames(data[[i]]) == "REF")
      alt.len = sapply(alt, length)
      num.alleles = max(alt.len)
      # Make this matrix of size num.alleles + 1 so that missing values can
      # be set to num.alleles + 1.
      repl = matrix("N", nrow = nrow(data[[i]]), ncol = num.alleles + 1)
      repl[,1] = data[[i]][,ref.col]
	  
      for(j in 1:num.alleles) {
        rng = which(j <= alt.len)
        repl[rng,j+1] = sapply(alt[rng], "[", j)
      } # for(j)
 
      # Replace missing values with "<max.value>/<max.value>"
      geno[geno == "."] = paste(num.alleles, num.alleles, sep = "/")
 
      # Convert the genotypes to a numeric matrix.
      geno = apply(geno, 2, function(z){ 
               matrix(unlist(strsplit(z, split = "/")), ncol = 2, byrow = TRUE)
             })
      g2 = matrix(0, nrow(qual), 2 * ncol(qual))
      g2[,2 * 1:ncol(geno) - 1] = as.numeric(geno[1:nrow(qual),])
      g2[,2 * 1:ncol(geno)] = as.numeric(geno[(nrow(qual)+1):nrow(geno),])
      geno = g2    
      g2 = apply(geno, 2, function(z) {
             repl[matrix(c(1:nrow(geno), z + 1), ncol = 2),drop = FALSE]
           })
		   
      # This is neccessary when there is only one SNP and R turns g2
      # into a vector.
      if(!is.matrix(g2)) {
        g2 = matrix(g2, nrow = 1)
      } # if(!is.matrix(g2))
	  
      # Create a two nucleotide character genotype matrix.
      g2 = split.data.frame(t(g2), factor(rep(1:ncol(qual), each = 2)))
      g2 = lapply(g2, function(z) { paste(z[1,], z[2,], sep = "") })
      geno = matrix(unlist(g2), nrow = nrow(qual), ncol = ncol(qual),
                    dimnames = dimnames(qual))
      rm(g2)
      gc()
	  
    } # else if(return.val == "allele")
	
    # Add the positions and alleles and combine the genotypes and quality
    # scores.
    pos = data[[i]][,1:5,drop = FALSE]
    info = data[[i]][,colnames(data[[i]]) == "INFO"]
	
    if(return.qual) {
	
      data[[i]] = data.frame(geno, qual)
      # This reorders the data with genotype and quality for each strain
      # next to each other.
      index = rep(1:(ncol(data[[i]])/2), each = 2)
      index[2 * 1:(ncol(data[[i]])/2)] = index[2 * 1:(ncol(data[[i]])/2)] + 
                                         (ncol(data[[i]]) / 2)
      data[[i]] = data[[i]][,index,drop = FALSE]
      colnames(data[[i]]) = sub("\\.1$", "qual", colnames(data[[i]]))
	  
    } else {
	
      data[[i]] = data.frame(geno, stringsAsFactors = FALSE)
	  
    } # else
	
    data[[i]] = cbind(pos, data[[i]])
	
    # If the user has requested variants consequences, add them in the last 
    # column.
    if(csq) {
      conseq = rep("", nrow(data[[i]]))
      which.csq = grep("CSQ", info)
      if(length(which.csq) > 0) {
        tmp = strsplit(info[which.csq], split = ";")
        tmp = sapply(tmp, function(z) { z[grep("^CSQ", z)] })
        conseq[which.csq] = tmp
        data[[i]] = data.frame(data[[i]], conseq = conseq)
      } # if(length(which.csq) > 0)
    } # if(csq)
    data[[i]]$POS = as.numeric(as.character(data[[i]]$POS))
  } # for(i)
  return(data)
  
} # read.vcf.snp()


################################################################################
# Retrieve Indels from a Sanger VCF file.
################################################################################
read.vcf.indel = function(data, keep, colnames, csq, return.val, return.qual) {

  data.not.empty = which(sapply(data, length) > 0)

  for(i in data.not.empty) {

    data[[i]] = strsplit(data[[i]], split = "\t")
    data[[i]] = matrix(unlist(data[[i]]), nrow = length(data[[i]]),
                byrow = TRUE, dimnames = list(NULL, colnames))
    data[[i]] = data[[i]][,keep,drop = FALSE]
    
    # Split up the genotype calls and quality scores.
    d2 = apply(data[[i]][,10:ncol(data[[i]]),drop=FALSE], 2, strsplit, split = ":")
    d2 = lapply(d2, function(z) { lapply(z, function(a) { a[c(1,6)] }) })
    geno = matrix(unlist(lapply(d2, function(z) { sapply(z, function(a) { a[1] }) })),
           ncol = length(d2), dimnames = list(NULL, names(d2)))
    qual = matrix(unlist(lapply(d2, function(z) { sapply(z, function(a) { a[2] }) })),
           ncol = length(d2), dimnames = list(NULL, names(d2)))
    qual[is.na(qual)] = "0"
    rm(d2)
    gc()
    # NOTE: MISSING VALUES AND REFERENCE ALLELE CALLS ARE BOTH DENOTED AS ".".
    # geno[geno == "."] = paste(num.alleles, num.alleles, sep = "/")
    # Replace "." calls, which may be reference or missing with "0/0".
    geno[geno == "."] = "0/0"
    
    if(return.val == "number") {
      geno = sub("/", "", geno)
      geno = sub("\\.", NA, geno)
      geno = matrix(as.numeric(geno), nrow(geno), ncol(geno), dimnames = 
             dimnames(geno))
    } else if(return.val == "allele") {
      # Create a matrix that maps numeric allele calls to nucleotide allele
      # calls.
      alt = strsplit(data[[i]][,colnames(data[[i]]) == "ALT"], split = ",")
      ref.col = which(colnames(data[[i]]) == "REF")
      alt.len = sapply(alt, length)
      num.alleles = max(alt.len)
      repl = matrix(NA, nrow = nrow(data[[i]]), ncol = num.alleles + 1)
      repl[,1] = data[[i]][,ref.col]
      for(j in 1:3) {
        rng = which(j <= alt.len)
        repl[rng,j+1] = sapply(alt[rng], function(z) { z[j] })
      } # for(j)

      # Convert the genotypes to a numeric matrix.
      geno = apply(geno, 2, function(z){ 
               matrix(unlist(strsplit(z, split = "/")), ncol = 2, byrow = TRUE)
             })
      g2 = matrix(0, nrow(qual), 2 * ncol(qual))
      g2[,2 * 1:ncol(geno) - 1] = as.numeric(geno[1:nrow(qual),])
      g2[,2 * 1:ncol(geno)] = as.numeric(geno[(nrow(qual)+1):nrow(geno),])
      geno = g2    
      g2 = apply(geno, 2, function(z) {
             repl[matrix(c(1:nrow(geno), z + 1), ncol = 2)]
           })

      # This is neccessary when there is only one SNP and R turns g2
      # into a vector.
      if(!is.matrix(g2)) {
        g2 = matrix(g2, nrow = 1)
      } # if(!is.matrix(g2))

      # Create a two allele character genotype matrix.
      # If all indels were homozygotes, we could just insert a single copy
      # of the allele. But there are het calls, so we concatenate the two
      # allele calls.
      g2 = split.data.frame(t(g2), factor(rep(1:ncol(qual), each = 2)))
      g2 = lapply(g2, function(z) { paste(z[1,], z[2,], sep = "") })
      geno = matrix(unlist(g2), nrow = nrow(qual), ncol = ncol(qual),
                    dimnames = dimnames(qual))
      rm(g2)
      gc()
    } # else if(return.val == "allele")

    # Add the positions and alleles and combine the genotypes and quality
    # scores.
    pos = data[[i]][,1:5]
    info = data[[i]][,colnames(data[[i]]) == "INFO"]
    data[[i]] = data.frame(geno, qual)
    index = rep(1:(ncol(data[[i]])/2), each = 2)
    index[2 * 1:(ncol(data[[i]])/2)] = index[2 * 1:(ncol(data[[i]])/2)] + (ncol(data[[i]]) / 2)
    data[[i]] = data[[i]][,index]
    colnames(data[[i]]) = sub("1$", "qual", colnames(data[[i]]))
    data[[i]] = cbind(pos, data[[i]])

    # If the user has requested variants consequences, add them in the last 
    # column.
    if(csq) {
      conseq = rep("", nrow(data[[i]]))
      which.csq = grep("CSQ", info)
      if(length(which.csq) > 0) {
        tmp = strsplit(info[which.csq], split = ";")
        tmp = sapply(tmp, function(z) { z[grep("^CSQ", z)] })
        conseq[which.csq] = tmp
        data[[i]] = data.frame(data[[i]], conseq = conseq)
      } # if(length(which.csq) > 0)
    } # if(csq)
  } # for(i)

  return(data)
  
} # read.vcf.indel()


################################################################################
# Retrieve structural variants from a Sanger VCF file.
################################################################################
read.vcf.sv = function(data, keep, colnames, csq, return.val, return.qual) {

  data.not.empty = which(sapply(data, length) > 0)

  for(i in data.not.empty) {

    data[[i]] = strsplit(data[[i]], split = "\t")
    data[[i]] = matrix(unlist(data[[i]]), nrow = length(data[[i]]),
                byrow = TRUE, dimnames = list(NULL, colnames))
    data[[i]] = data[[i]][,keep,drop = FALSE]
    
    # Split up the genotype calls and quality scores.
    d2 = apply(data[[i]][,5:ncol(data[[i]]),drop=FALSE], 2, strsplit, split = ";")
    d2 = lapply(d2, function(z) { lapply(z, function(a) { a[c(1,6)] }) })
    geno = matrix(unlist(lapply(d2, function(z) { sapply(z, function(a) { a[1] }) })),
           ncol = length(d2), dimnames = list(NULL, names(d2)))
    qual = matrix(unlist(lapply(d2, function(z) { sapply(z, function(a) { a[2] }) })),
           ncol = length(d2), dimnames = list(NULL, names(d2)))
    rm(d2)
    gc()

    if(return.val == "number") {

      geno = sub("/", "", geno)
      geno = sub("\\.", NA, geno)
      geno = matrix(as.numeric(geno), nrow(geno), ncol(geno), dimnames = 
             dimnames(geno))

    } else if(return.val == "allele") {

      # Create a matrix that maps numeric allele calls to nucleotide allele
      # calls.
      alt = strsplit(data[[i]][,colnames(data[[i]]) == "ALT"], split = ",")
      ref.col = which(colnames(data[[i]]) == "REF")
      alt.len = sapply(alt, length)
      num.alleles = max(alt.len)
      # Make this matrix of size num.alleles + 1 so that missing values can
      # be set to num.alleles + 1.
      repl = matrix(NA, nrow = nrow(data[[i]]), ncol = num.alleles + 1)
      repl[,1] = data[[i]][,ref.col]
      for(j in 1:num.alleles) {
        rng = which(j <= alt.len)
        repl[rng,j+1] = sapply(alt[rng], function(z) { z[j] })
      } # for(j)
 
      # Replace missing values with "<max.value>/<max.value>"
      geno[geno == "."] = paste(num.alleles, num.alleles, sep = "/")
 
      # Convert the genotypes to a numeric matrix.
      geno = apply(geno, 2, function(z){ 
               matrix(unlist(strsplit(z, split = "/")), ncol = 2, byrow = TRUE)
             })
      g2 = matrix(0, nrow(qual), 2 * ncol(qual))
      g2[,2 * 1:ncol(geno) - 1] = as.numeric(geno[1:nrow(qual),])
      g2[,2 * 1:ncol(geno)] = as.numeric(geno[(nrow(qual)+1):nrow(geno),])
      geno = g2    
      g2 = apply(geno, 2, function(z) {
             repl[matrix(c(1:nrow(geno), z + 1), ncol = 2)]
           })
      # This is neccessary when there is only one SNP and R turns g2
      # into a vector.
      if(!is.matrix(g2)) {
        g2 = matrix(g2, nrow = 1)
      } # if(!is.matrix(g2))
      # Create a two nucleotide character genotype matrix.
      g2 = split.data.frame(t(g2), factor(rep(1:ncol(qual), each = 2)))
      g2 = lapply(g2, function(z) { paste(z[1,], z[2,], sep = "") })
      geno = matrix(unlist(g2), nrow = nrow(qual), ncol = ncol(qual),
                    dimnames = dimnames(qual))
      rm(g2)
      gc()
    } # else if(return.val == "allele")

    # Add the positions and alleles and combine the genotypes and quality
    # scores.
    pos = data[[i]][,1:5]
    info = data[[i]][,colnames(data[[i]]) == "INFO"]
    data[[i]] = data.frame(geno, qual)
    index = rep(1:(ncol(data[[i]])/2), each = 2)
    index[2 * 1:(ncol(data[[i]])/2)] = index[2 * 1:(ncol(data[[i]])/2)] + 
                                       (ncol(data[[i]]) / 2)
    data[[i]] = data[[i]][,index]
    colnames(data[[i]]) = sub("1$", "qual", colnames(data[[i]]))
    data[[i]] = cbind(pos, data[[i]])

    # If the user has requested variants consequences, add them in the last 
    # column.
    if(csq) {
      conseq = rep("", nrow(data[[i]]))
      which.csq = grep("CSQ", info)
      if(length(which.csq) > 0) {
        tmp = strsplit(info[which.csq], split = ";")
        tmp = sapply(tmp, function(z) { z[grep("^CSQ", z)] })
        conseq[which.csq] = tmp
        data[[i]] = data.frame(data[[i]], conseq = conseq, stringsAsFactors = FALSE)
      } # if(length(which.csq) > 0)
    } # if(csq)

  } # for(i)

  return(data)
  
} # read.vcf.sv()



################################################################################
# Get the strain names available in the requested VCF file.
################################################################################
get.vcf.strains = function(vcf.file) {
  # Query Tabix indexed VCF file.
  tabix = TabixFile(file = vcf.file)
  open(con = tabix)
  hdr = headerTabix(file = tabix)
  close(con = tabix)
  # Split up the strains.
  strains = sub("^#", "" ,hdr$header[length(hdr$header)])
  strains = strsplit(strains, split = "\t")[[1]]
  if(length(grep("snp", vcf.file)) > 0) {
    retval = strains[-1:-9]
  } else if(length(grep("indel", vcf.file)) > 0) {
    retval = strains[-1:-9]
  } else if(length(grep("SV", vcf.file)) > 0) {
    retval = strains[-1:-4]
  } else {
    stop(paste("Unknown file type. The VCF file name must contain one of",
         "snp, indel or SV in order to identify the type of file to parse."))
  } # else
  return(retval)
} # get.vcf.strains()
