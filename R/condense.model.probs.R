################################################################################
# Given the HMM 36 state genotype probability files, condense them to 8 founder
# states and optionally write out the large 3D array to a *.Rdata file.
# This requires that the genotype.probs.txt files have been converted to 
# genotype.probs.Rdata files.
# Daniel Gatti
# Dan.Gatti@jax.org
# Jan. 31, 2012
# Aug. 4, 2013: Added dominance and full model arrays.
################################################################################
# Arguments: path: character, full path to the directory where the
#                  *.genotype.probs.Rdata files are.
#            write: character, if missing, writes nothing and returns
#                   the large 3D array. Otherwise, the full path to a file to 
#                   write the founder probabilities out to.
#            model: character, one of "additive", "dominance" or "full", 
#                   indicating the type of array to return. 'Additive' condenses
#                   the 36 genotypes down to 8 founder allele dosages. 
#                   'Dominance' condenses them down to 16 values, the first 8
#                   are additive values and the last 8 are dominance values.
#                   'Full' returns all 36 states.
#            cross: character containing the cross type. Typically "DO, "CC, 
#                   "HS" or "DOF1".
condense.model.probs = function(path = ".", write = "founder.probs.Rdata",
                       model = c("additive", "dominance", "full"),
                       cross = "DO") {

  model = match.arg(model)
  path  = add.slash(path)

  if(!is.character(path)) {
    stop(paste("\\'path\\' must be a characer vector containing the path to",
        "the directory where the genoprobs files are."))
  } # if(!is.character(path))

  # Get the files that we need to process.
  files = dir(path, pattern = ".genotype.probs.Rdata$", full.names = TRUE)

  # Get their sample IDs.
  samples = dir(path, pattern = "genotype.probs.Rdata$", full.names = FALSE)
  samples = sub("\\.genotype\\.probs\\.Rdata$", "", samples)
  model.probs = NULL

  if(model == "additive") {
    model.probs = get.additive(files = files, samples = samples)
  } else if(model == "dominance") {
    model.probs = get.dominance(files = files, samples = samples)
  } else if(model == "full") {
    model.probs = get.full(files = files, samples = samples)
  } else {
    stop(paste("condense.model.probs:Invalid model '", model,"'. Please",
         "enter one of additive, dominance or full in the model argument."))
  } # else

  attr(model.probs, "cross") = cross

  print("Writing data...")
  save(model.probs, file = write)

} # condense.model.probs()


# Helper function to get the additive model probabilities, which condenses
# the 36 genotypes into additive allele dosages.
get.additive = function(files, samples) {

  # Create a large 3D array with samples/states/SNPs in dimensions 1,2,3.
  print(files[1])
  load(files[1]) # load in prsmth
  prsmth = as.matrix(prsmth)
  founders = sort(unique(unlist(strsplit(colnames(prsmth), split = ""))))
  model.probs = array(0, c(length(files), length(founders), nrow(prsmth)),
                dimnames = list(samples, founders, rownames(prsmth)))

  # Create a matrix that we can use to multiply the 36 state probabilities.
  mat = get.diplotype2haplotype.matrix(colnames(prsmth))

  # Place the founder probabilities in the matrix.
  model.probs[1,,] = t(prsmth %*% mat)
  for(i in 2:length(files)) {
    print(files[i])
    load(files[i]) # load in prsmth
    model.probs[i,,] = t(prsmth %*% mat)
  } # for(i)

  model.probs

} # get.additive()


# Helper function to get the additive model probabilities, which condenses
# the 36 genotypes into additive and dominant dosages.
get.dominance = function(files, samples) {
  # Create a large 3D array with samples/states/SNPs in dimensions 1,2,3.
  load(files[1]) # load in prsmth
  prsmth = t(prsmth)
  founders = sort(unique(unlist(strsplit(rownames(prsmth), split = ""))))
  model.probs = array(0, c(length(files), 2 * length(founders), ncol(prsmth)),
                dimnames = list(samples, paste(rep(founders, 2), rep(c("add", "dom"), 
                each = length(founders), sep = ".")), colnames(prsmth)))
  # Create a matrix that we can use to multiply the 36 state probabilities
  # to get additive values.
  addmat = matrix(0, length(founders), nrow(prsmth), dimnames = list(
               founders, rownames(prsmth)))
  m = sapply(founders, grep, colnames(addmat))
  for(i in 1:nrow(addmat)) {
    addmat[i,m[i,]] = 0.5
  } # for(i)
  addmat[matrix(c(1:length(founders), diag(m)), ncol = 2)] = 1
  # Create a matrix that we can use to multiply the 36 state probabilities
  # to get additive values.
  dommat = matrix(0, length(founders), nrow(prsmth), dimnames = list(
               founders, rownames(prsmth)))
  m = sapply(founders, grep, colnames(dommat))
  for(i in 1:nrow(dommat)) {
    dommat[i,m[i,]] = 1.0
  } # for(i)
  dommat[matrix(c(1:length(founders), diag(m)), ncol = 2)] = 1
  # Place the founder probabilities in the matrix.
  model.probs[1,1:length(founders),] = addmat %*% prsmth
  model.probs[1,(length(founders) + 1):dim(model.probs)[[2]],] = dommat %*% prsmth
  for(i in 2:length(files)) {
    print(files[i])
    load(files[i]) # loadin prsmth
    prsmth = t(prsmth)
    model.probs[i,1:length(founders),] = addmat %*% prsmth
    model.probs[i,(length(founders) + 1):dim(model.probs)[[2]],] = dommat %*% prsmth
  } # for(i)
 
  model.probs

} # get.dominance()


# Helper function to get the full model probabilities, which just gathers
# the 36 states for each samples into one file.
get.full = function(files, samples) {

  # Create a large 3D array with samples/states/SNPs in dimensions 1,2,3.
  load(files[1]) # load in prsmth
  prsmth = t(prsmth)
  model.probs = array(0, c(length(files), nrow(prsmth), ncol(prsmth)),
                dimnames = list(samples, rownames(prsmth), colnames(prsmth)))

  # Place the probabilities in the matrix.
  model.probs[1,,] = prsmth

  for(i in 2:length(files)) {
    print(files[i])
    load(files[i]) # loadin prsmth
    prsmth = t(prsmth)
    model.probs[i,,] = prsmth
  } # for(i)

  model.probs

} # get.full()


get.diplotype2haplotype.matrix = function(states) {

  spl = strsplit(states, split = "")
  names(spl) = states
  spl = lapply(spl, factor, levels = sort(unique(unlist(spl))))
  spl = lapply(spl, table)
  mat = matrix(unlist(spl), length(spl[[1]]), length(spl), dimnames = 
        list(names(spl[[1]]), names(spl)))
  mat = t(mat * 0.5)

} # get.diplotype2haplotype.matrix()

