################################################################################
# Verify that the sample matrices have data in them and that they are the same
# dimensions.
# Arguments: x: numeric matrix of normalized X intensities with samples in rows
#               and SNPs in columns. Both dimensions must be named.
#            y: numeric matrix of normalized Y intensities with samples in rows
#               and SNPs in columns. Both dimensions must be named.
#            is.founder.F1: boolean vector that is TRUE if a sample is a 
#               founder or F1. Length must equal number of rows in theta & rho.
#            init.means: 3D array of initial theta and rho mean values to use.
#                        (num states x num SNPs x 2). May be NULL if founders &
#                        F1s are provided with the data. 
#            init.covars: 3D array of initial theta and rho covariance values to
#                         use. (num states x num SNPs x 2). May be NULL if
#                         founders & F1s are provided with the data.
#            chr.to.run: character vector with the Chr to run the HMM on.
#            snps: data.frame with the SNPs to be used in the analysis.
#            sex: character vector with sex IDs which we need for the X Chr.
#            plot: boolean that is true if the user would like to plot a sample
#                  chromosome as the model progresses.
check.data.for.errors = function(x, y, is.founder.F1, init.means,
                                 init.covars, chr.to.run, snps, sex, plot) {
  # Null X data.
  if(is.null(x)) {
    stop("x is null. Please enter a non-null x matrix.")
  } # if(is.null(x))

  # Null Y data.
  if(is.null(y)) {
    stop("y is null. Please enter a non-null y matrix.")
  } # if(is.null(y))

  # Null founder data.
  if(is.null(is.founder.F1)) {
    stop("is.founder.F1 is null. Please enter non-null founder & F1 data.")
  } # if(is.null(is.founder.F1))

  # Founder IDs are all two letters long.
  if(!all(nchar(rownames(x)[is.founder.F1]) == 2)) {
    stop(paste("All of the founder IDs are not two letters long. They must",
         "all be two letter codes. The problem names are:", paste(
         rownames(x)[is.founder.F1][nchar(rownames(x)[is.founder.F1]) != 2],
         collapse = ",")))
  } # if(!all(nchar(rownames(x)[is.founder.F1])))
   
  # X & Y dimensions equal.
  if(any(dim(x) != dim(y))) {
    stop(paste("The dimensions of the x (", nrow(x), ", ", ncol(x) ,
         ") and y (", nrow(y), ", ", ncol(y),
         ") matrices differ. Please enter x and y matrices with the",
         " same dimensions.", sep = ""))
  } # if(any(dim(x) != dim(y))

  # Rownames and colnames match in X & Y.
  if(!all(rownames(x) == rownames(y))) {
    stop(paste("The rownames(x) and rownames(y) are not equal. Please verify",
         "that both X & Y have the same rownames."))
  } # if(!all(rownames(x) == rownames(y)))

  if(!all(colnames(x) == colnames(y))) {
    stop(paste("The colnames(x) and colnames(y) are not equal. Please verify",
         "that both X & Y have the same colnames."))
  } # if(!all(colnames(x) == colnames(y)))

  # Founder IDs == num.samples 
  if(length(is.founder.F1) != nrow(x)) {
    stop(paste("The number of samples in the founder IDs (",
         length(is.founder.F1),
         ") does not match the number of samples in the data (", nrow(x),
         ").", sep = ""))
  } # if(nrow(is.founder.F1) != nrow(x))

  # Some founder IDs true or init.means non-NULL with the correct dimensions.
  if(sum(is.founder.F1) == 0) {
  	if(is.null(init.means)) {
      stop(paste("None of the founder IDs are TRUE.",
          "Please include founder and F1 data in the samples or include state",
          "mean and varinace estimates in init.means and init.covars."))
    } else {
      # Verify that init.means has the correct dimensions.
      if(ncol(init.means) != ncol(x)) {
        stop(paste("The init.means and X matrices have differing numbers of",
              "SNPs. Please make sure that they have the same dimensions."))
      } # if(ncol(init.means) != ncol(x)) 
      
     if(ncol(init.covars) != ncol(x)) {
        stop(paste("The init.covars and X matrices have differing numbers of",
              "SNPs. Please make sure that they have the same dimensions."))
      } # if(ncol(init.covars) != ncol(x)) 
    } # else
  } # if(sum(is.founder.F1) == 0)  

  # No NA values in X.
  if(any(is.na(x))) {
    print(paste("There are ", sum(is.na(x)), " NA values in the X matrix.",
          "Please impute or replace these with floating point numbers."))
    return(NULL)
  } # if(any(is.na(x)))

  # No NA values in Y.
  if(any(is.na(y))) {
    print(paste("There are ", sum(is.na(y)), " NA values in the Y",
          "matrix. Please impute or replace these with floating point numbers."))
    return(NULL)
  } # if(any(is.na(y)))

  # If the user has specified any sex chromosomes to run (X, Y, M), then verify
  # that we have sex information.
  if(any(c("X", "Y", "M") %in% chr.to.run)) {
  	if(is.null(sex)) {
      print(paste("In order to analyze the sex chromosomes, the sex variable",
            "must have M & F values for each sample."))
      return(NULL)
  	} # if(sex == NULL)
  } # if(any(c("X", "Y", "M") %in% chr.to.run))

  # Verify that, if the X chr is being run, we have examples of the all of the
  # founders with the correct sex.
  sex = toupper(sex)
  if("X" %in% chr.to.run) {
    founders = sort(unique(unlist(strsplit(rownames(x)[is.founder.F1], split = ""))))
    founders = paste(founders, founders, sep = "")
    if(any(sex == "F")) {
      female.founders = is.founder.F1 & rownames(x) %in% founders & sex == "F"
      female.founders = rownames(x)[female.founders]
      if(!all(founders %in% female.founders)) {
        stop(paste("If you are genotyping the X chromosome with female samples,",
             "you must provide at least one female sample from each founder strain."))
      } # if(!all(founders %in% female.founders))
    } #if(any(sex == "F"))

    if(any(sex == "M")) {
      male.founders = is.founder.F1 & rownames(x) %in% founders & sex == "M"
      male.founders = rownames(x)[male.founders]
      if(!all(founders %in% male.founders)) {
        stop(paste("If you are genotyping the X chromosome with male samples,",
             "you must provide at least one male sample from each founder strain."))
      } # if(!all(founders %in% male.founders))
    } # if(any(sex == "M"))
  } # if("X" %in% chr.to.run)

  # Provide a warning if the sex of any samples is unclear.
  predicted.sex = predict.sex(x, y, snps, F)
  sex.mismatch = which(sex != predicted.sex)
  if(any(sex.mismatch)) {
    warning(paste("The following samples may not have the correct sex",
            "assignment:", paste(names(predicted.sex)[sex.mismatch],
            collapse = ",")))
  } # if(any(sex.mismatch))

  # Verify that the SNP IDs line up everywhere.
  if(!all(snps[,1] == colnames(x))) {
    stop(paste("The SNP IDs in X & Y must match the SNP IDs in the package."))
  } # if(!all(snps[,1] == colnames(x)))
  
  return("OK")
} # check.data.for.errors()

