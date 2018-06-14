# Given model parameters, run forward-backward and output genoprobs.
# data: list containing elements named x, y, sex and gen. x & y are numeric matrices
#       with num.samples rows and num.markers columns that contain X and Y intensities.
#       sex is a named character vector containing "F" or "M". gen is a numeric vector
#       containing the DO outbreeding generation.
# param.files: list containing the full path to the parameter files. These have the form 
#              *.final.means.Rdata and *.final.covars.Rdata. Each list element should be
#              a character vector with two file names in it. Each list must be names
#              with the chromosome or X.male and X.female.
# markers: data.frame containing the marker information. Marker ID, Chr, Mb postion and
#          cM in columns 1 through 4.
# NOTE: For the paralell version to work, you must have openmpi, Rmpi, snow and 
#       doParallel installed. Also, you must qsub with the number of nodes = to the 
#       number of chromosomes.
calc.genoprob2 = function(data, param.files, markers,
                 cross = c("DO", "CC", "HS", "DOF1", "other"),
                 output.dir = ".") {
  if(missing(data)) {
    stop("Required argument \'data\' argument missing.")
  } # if(missing(data))
  if(missing(param.files)) {
    stop("Required argument \'param.files\' argument missing.")
  } # if(missing(data))
  if(missing(markers)) {
    stop("Required argument \'markers\' argument missing.")
  } # if(missing(markers))
  # Split the markers by chromosome.
  markers = split(markers, markers[,2])
  markers = markers[order(as.numeric(names(markers)))]
  cross = match.arg(cross)
  output.dir = add.slash(output.dir)
  # Sort the param files by chromosome.
  if(is.null(names(param.files))) {
    stop(paste("param.files must have the chromosome names in names(param.files)."))
  } else {
    param.files = param.files[order(as.numeric(names(param.files)))]
  } # else
  # Split the data up by chromosome. Separate males and females on X. We will place 
  # the file name for the parameter files in data and let the clusters load in the data.
  newdata = vector("list", length(markers))
  names(newdata) = names(markers)
  # Split the X chromosome into females and males.
  if(any(names(markers) == "X")) {
    if(any(data$sex == "F")) {
      names(newdata)[[length(newdata)]] = "X.female"
      if(any(data$sex == "M")) {
        newdata = c(newdata, "")
        names(newdata)[[length(newdata)]] = "X.male"
      } # 
    } else {
      names(newdata)[[length(newdata)]] = "X.male"
    } # else
  } # if(any(names(markers) == "X"))
  auto = which(!is.na(as.numeric(names(markers))))
  if(length(auto) > 0) {
    for(chr in auto) {
      newdata[[chr]] = list(x = data$x[,markers[[chr]][,1]],
                            y = data$y[,markers[[chr]][,1]],
                            sex = data$sex,
                            gen = data$gen,
                            means  = param.files[[chr]][1],
                            covars = param.files[[chr]][2],
                            markers = markers[[chr]],
                           out.file = paste0(output.dir, "chr", chr, "_genoprobs.Rdata"))
    } # for(i)
  } # if(length(auto) > 0)
  # X chromosome.
  females = which(toupper(data$sex) == "F")
  if(length(females) > 0) {
    newdata[["X.female"]] = list(x = data$x[females, markers[["X"]][,1]],
                                 y = data$y[females, markers[["X"]][,1]],
                                 sex = data$sex[females],
                                 gen = data$gen[females],
                                 means  = param.files[["X.female"]][1],
                                 covars = param.files[["X.female"]][2],
                                 markers = markers[["X"]],
                                 out.file = paste0(output.dir, "chrX_female_genoprobs.Rdata"))
  } #if(length(females) > 0)
  males = which(toupper(data$sex) == "M")
  if(length(males) > 0) {
    newdata[["X.male"]] = list(x = data$x[males, markers[["X"]][,1]],
                               y = data$y[males, markers[["X"]][,1]],
                               sex = data$sex[males],
                               gen = data$gen[males],
                               means  = param.files[["X.male"]][1],
                               covars = param.files[["X.male"]][2],
                               markers = markers[["X"]],
                               out.file = paste0(output.dir, "chrX_male_genoprobs.Rdata"))
  } # if(length(males) > 0)
  # Remove unused data.
  rm(data, param.files, markers, males, females)
  gc()
  # Create the cluster. We default to one cluster per chromosome (plus one if we have
  # the X chromosome and two sexes).
  cl = makeCluster(spec = 5, outfile = "cluster_out.txt")
  registerDoParallel(cl = cl)
  clusterEvalQ(cl = cl, expr = library(DOQTL))
  # For each chromosome, calculate the genoprobs on the autosomes and 
  # write out the results as 3-dimensional arrays in *.Rdata files.
  result = foreach(i = iter(newdata)) %dopar% {
    DOQTL:::genoprob.helper(i)
  } # for(i)
  samples = result[[1]][[2]]
  states  = result[[1]][[1]]
  # We need to read in the male and female probs and combine them.
  if(length(grep("X", names(newdata))) > 0) {
    message("Merging X chromosome results ...")
    fem.probs = NULL
    mal.probs = NULL
    x.markers = NULL
    if(any(names(newdata) == "X.female")) {
      load(newdata[["X.female"]]$out.file)
      fem.probs = probs
      x.markers = dimnames(fem.probs)[[3]]
    } # if(any(names(newdata) == "X.female"))
    if(any(names(newdata) == "X.male")) {
      load(newdata[["X.male"]]$out.file)
      mal.probs = probs
      rownames(mal.probs) = paste0(rownames(mal.probs), rownames(mal.probs))
      x.markers = dimnames(mal.probs)[[3]]
    } # if(any(names(newdata) == "X.male"))
    # Merge the X chromosome results.
    probs = array(0, c(length(states), 
                  length(samples), length(x.markers)), dimnames = 
                  list(states, samples, x.markers))
    if(any(names(newdata) == "X.female")) {
      probs[,match(colnames(fem.probs), colnames(probs)),] = fem.probs
      rm(fem.probs)
    } # if(any(names(newdata) == "X.female"))
    if(any(names(newdata) == "X.male")) {
      probs[match(rownames(mal.probs), rownames(probs)),
            match(colnames(mal.probs), colnames(probs)),] = mal.probs
      rm(mal.probs)
    } # if(any(names(newdata) == "X.male"))
    save(probs, file = paste0(output.dir, "chrX_genoprobs.Rdata"))
  } # if(length(grep("X", names(newdata))) > 0)
  stopCluster(cl)
} # calc.genoprob2()
# Function to calculate the genoprobs for each chromosome.
# The data argument is a list with eight named elements:
#  x: numeric matrix containing the X intensities. Samples in rows, markers 
#     in columns. Rownames and colnames must contain sample IDs and mrkers IDs.
#  y: numeric matrix containing the Y intensities. Samples in rows, markers 
#     in columns. Rownames and colnames must contain sample IDs and mrkers IDs.
#  sex: character vector containing M or F to indicate sex. Samples IDs must 
#       be in names.
#  gen: numeric vector containing the DO outbreeding generation. Samples IDs must 
#       be in names.
#  means: character string containing the full path to the cluster means file
#         for the current chromosome.
#  covars: character string containing the full path to the cluster covariance file
#         for the current chromosome.
#  markers: data.frame containing at least 4 columns with marker ID, chr, Mb postion
#           and cM postion, respectively.
# out.file: character string with the full path of the *.Rdata file to write out.
genoprob.helper = function(data, out.file) {
  # Load in the paramaters.
  load(data$means)
  data$means = theta.rho.means
  load(data$covars)
  data$covars = theta.rho.covars
  # Synchronize the samples.
  keep = intersect(rownames(data$x), names(data$sex))
  message(paste("Chr", data$markers[1,2], ": using the", length(keep), 
          "samples that intersect between x, y, sex and gen."))
  data$x = data$x[keep,,drop = FALSE]
  data$y = data$y[keep,,drop = FALSE]
  data$sex = data$sex[keep]
  data$gen = data$gen[keep]
  # Synchronize the markers.
  keep = intersect(intersect(colnames(data$x), dimnames(data$means)[[3]]), data$markers[,1])
  message(paste("Chr", data$markers[1,2], ": using the", length(keep), 
          "markers that intersect between x, y, means, covars and markers."))
  data$x = data$x[,keep,drop = FALSE]
  data$y = data$y[,keep,drop = FALSE]
  data$means  = data$means[,,keep, drop = FALSE]
  data$covars = data$covars[,,keep, drop = FALSE]
  data$markers = data$markers[data$markers[,1] %in% keep,]
  # Convert X and Y to theta and rho.
  data = list(theta = (2.0 / pi) * atan2(data$y, data$x), 
              rho = sqrt(data$x^2 + data$y^2), 
              sex = data$sex, 
              gen = data$gen,
              means = data$means,
              covars = data$covars,
              markers = data$markers,
              out.file = data$out.file)
  # Get transition probabilities.
  # For DO samples, we have different transition probabilities for each
  # generation.
  a = do.trans.probs(states = rownames(data$means), snps = data$markers, 
                     chr = data$markers[1,2], sex = data$sex[1], gen = data$gen)
  # Get emission probabilities for each sample.
  b = emission.probs.intensity(data = data, params = list(r.t.means = data$means, 
      r.t.covars = data$covars))
  # Get initial probs and allocate memory.
  init.hmm = initialize.hmm(snps = colnames(data$rho), samples = 
             rownames(data$rho), states = rownames(data$means))
  # Run filtering and smoothing.
  for(i in 1:length(a)) {
    # Run foreward/backward algorithm.
    print(paste("gen", names(a)[i]))
    curr.gen = which(data$gen == names(a)[i])
    lltmp = 0
    res = .C(C_filter_smooth_intensity, 
             dims = as.integer(dim(init.hmm$prsmth[,curr.gen,,drop = FALSE])),
             a = as.double(a[[i]]), 
             b = as.double(b[,curr.gen,,drop = FALSE]),
             prsmth = as.double(init.hmm$prsmth[,curr.gen,,drop = FALSE]),
             init = as.double(init.hmm$init),
             loglik = as.double(lltmp))
    init.hmm$prsmth[,curr.gen,] = res$prsmth
  } # for(i)
  message("Writing Results...")
  probs = exp(init.hmm$prsmth)
  save(probs, file = data$out.file)
  # Return the genotype states, samples and markers.
  return(dimnames(probs))
} # genoprob.helper()
