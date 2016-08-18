################################################################################
# Perform haplogroup assignment on chr M or Chr Y in DO samples genotyped on 
# the MegaMUGA or GigaMUGA.
# Daniel Gatti
# dan.gatti@jax.org
# Aug. 16, 2016
################################################################################
# Arguments: x: numeric matrix containing X intensities. Samples in rows and 
#               markers in columns. Must be same dimension as y.
#            y: numeric matrix containing Y intensities. Samples in rows and 
#               markers in columns. Must be same dimension as x.
#         geno: character matrix containing the allele calls. Samples in rows
#               and markers in columns.
#        array: character string that is either "megamuga" or "gigamuga", 
#               indicating the array platform of the data. Default = "gigamuga".
#       method: character string that is either "intensity" or "allele",
#               indicating whether to use intensites or allele calls to perform
#               the haplogroup calls.
# If method = "intensity", only x and y are required.  If method = "allele",
# only geno is required.
# We retrieve the founder data from ftp://ftp.jax.org/MUGA. 
get.haplogroup = function(x, y, geno, chr = c("M", "Y"), array = c("gigamuga",
                 "megamuga"), method = c("intensity", "allele")) {

  array = match.arg(array)
  chr = match.arg(chr)
  method = match.arg(method)
  
  # Verify that we have intensites for the intensity method or genotypes for
  # the allele method.
  if(method == "intensity") {

    if(missing(x) | missing(y)) {
      stop("You must provide x and y intensites for the intensity method.")
    } # if(missing(x) | missing(y))

    founders = read.muga.data(array = array, method = method, sampletype = "DO")
    return(get.haplogroup.intensity(x = x, y = y, founders = founders, 
           chr = chr, array = array))

  } else if(method == "allele") {

    if(missing(geno)) {
      stop("You must provide geno for the allele method.")
    } # if(missing(geno))

    founders = read.muga.data(array = array, method = method, sampletype = "DO")
    return(get.haplogroup.allele(geno = geno, chr = chr, array = array))

  } # else

  stop("get.haplogroup: method was neither intensity or allele.")

} # get.haplogroup()


get.haplogroup.intensity = function(x, y, founders, chr, array) {

  snps = founders$snps
  founders = founders[names(founders) != "snps"]

  # Return value.
  cl = NULL

  if(chr == "Y") {

    # Get the sex of the samples and keep only the males.
    # Keep only the Chr Y markers that have well-separated clusters 
    # (snps$haploChrY).
    sex = sex.predict(x = x, y = y, snps = snps, plot = FALSE)
    males = which(sex == "M")
    x = x[males, snps$haploChrY]
    y = y[males, snps$haploChrY]
    founders = keep.homozygotes(founders)
    fx = founders$x[founders$sex == "M", snps$haploChrY]
    fy = founders$y[founders$sex == "M", snps$haploChrY]
    rownames(fx) = rownames(fy) = founders$code[founders$sex == "M"]
    # Run PCA on teh combined X and Y data.
    pcx = do.pca(x = rbind(x, fx), y = rbind(y, fy))

    # Cluster using PAMK with 6 clusters. We can identify these in the DO.
    cl = pamk(data = pcx$scores, krange = 6, usepam = FALSE)
    cl = cl$pamobject$clustering
    # Match the cluster numbers to the founders.
    founder.cl = cl[names(cl) %in% founders$code]
    cl = cl[!names(cl) %in% founders$code]
    founder.cl = split(founder.cl, founder.cl)
    founder.cl = lapply(founder.cl, function(z) { 
                   paste0(unique(names(z)), collapse = "") 
                 })
    # NOTE: we will call the CCBBEE cluster "BB".
    founder.cl = unlist(founder.cl)
    founder.cl[grep("BB", founder.cl)] = "BB"
    cl = founder.cl[cl]
    names(cl) = rownames(x)

  } else if(chr == "M") {

    # Keep only the Chr M markers that have well-separated clusters 
    # (snps$haploChrM).
    x = x[,snps$haploChrM]
    y = y[,snps$haploChrM]
    founders = keep.homozygotes(founders)
    fx = founders$x[,snps$haploChrM]
    fy = founders$y[,snps$haploChrM]
    rownames(fx) = rownames(fy) = founders$code
    # Run PCA on the combined X and Y data.
    pcx = do.pca(x = rbind(x, fx), y = rbind(y, fy))

    # Cluster using PAMK with 5 clusters. We can identify these in the DO.
    cl = pamk(data = pcx$scores, krange = 5, usepam = FALSE)
    cl = cl$pamobject$clustering
    # Match the cluster numbers to the founders.
    founder.cl = cl[names(cl) %in% founders$code]
    cl = cl[!names(cl) %in% founders$code]
    founder.cl = split(founder.cl, founder.cl)
    founder.cl = lapply(founder.cl, function(z) { 
                   paste0(unique(names(z)), collapse = "") 
                 })
    # NOTE: we will call the CCBBEE cluster "BB".
    founder.cl = unlist(founder.cl)
    founder.cl[grep("BB", founder.cl)] = "BB"
    cl = founder.cl[cl]
    names(cl) = rownames(x)

  } else {
    stop("get.haplogroup.intensity: chr was neither M nor Y.")
  } # else

  return(cl)

} # get.haplogroup.intensity()


get.haplogroup.allele = function(geno, chr, array) {

} # get.haplogroup.allele()


do.pca = function(x, y) {

  # Impute missing data.
  if(any(is.nan(range(x)))) {
    x = impute.knn(data = x)$data
  } # if(any(is.nan(range(x))))

  if(any(is.nan(range(y)))) {
    y = impute.knn(data = y)$data
  } # if(any(is.nan(range(y))))

  xy = cbind(x, y)
  xy = apply(xy, 2, scale)
  rownames(xy) = rownames(x)

  return(princomp(xy))

} # do.pca()

get.founder.probs = function(geno) {

  founders = unique(rownames(geno))
  retval = matrix(NA, length(founders), ncol(geno), dimnames =
           list())
  for(i in 1:length(founders)) {

  } # for(i)

} # get.founder.probs()


plot.haplogroups = function(pcx, cl) {

  col = do.colors
  col[1,3] = "#FFC000"
  rownames(col) = paste0(col[,1], col[,1])
  col = col[cl,3]
  layout(matrix(1:3, 1, 3))
  plot(pcx$scores[,c(2,1)], pch = 16, col = col)
  plot(pcx$scores[,c(3,1)], pch = 16, col = col)
  plot(pcx$scores[,c(4,1)], pch = 16, col = col)

} # plot.haplogroups()

