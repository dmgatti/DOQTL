################################################################################
# Predict the sex of the samples using the mean X and Y Chr intensity values for
# each sample.
# Arguments: x: numeric matrix, num.samples x num.snps, with X intensities for
#               all samples.
#            y: numeric matrix, num.samples x num.snps, with Y intensities for
#               all samples.
#            snps: data.frame with three columns: SNP ID, Chr, location.
#            plot: Boolean that will create a plot of mean X and Y Chr 
#                  intensities if TRUE.  Default = FALSE.
# Returns: Character vector with sex assignments based on linear discriminant
#          analysis.
sex.predict = function(x, y, snps, plot = FALSE) {
  if(all(snps[,2] != "X")) {
    stop(paste("There are no X chromosome SNPs in snps. X and Y chromosome",
        "SNPs are required to predict sex."))
  }
  if(all(snps[,2] != "Y")) {
    stop(paste("There are no Y chromosome SNPs in snps. X and Y chromosome",
        "SNPs are required to predict sex."))
  }
  # Convert X and Y to matrices, if needed.
  if(!is.matrix(x)) {
    x = as.matrix(x)
  } # if(!is.matrix(x))
  if(!is.matrix(y)) {
    y = as.matrix(y)
  } # if(!is.matrix(y))
  # Synch up the X & Y SNPs.
  x = x[,colnames(x) %in% snps[,1]]
  y = y[,colnames(y) %in% snps[,1]]
  snps = snps[match(colnames(x), snps[,1]),]
  # Keep the X and Y chromosome data.
  x.rng = which(snps[,2] == "X")
  x.int = rowMeans(x[,x.rng] + y[,x.rng], na.rm = TRUE)
  y.rng = which(snps[,2] == "Y")
  y.int = rowMeans(x[,y.rng] + y[,y.rng], na.rm = TRUE)
  keep = which(!(is.na(x.int) & is.na(y.int)))
  if(length(keep) < length(x.int)) {
    warning(paste("Removing", length(x.int) - length(keep),
            "samples with NaN intensity values."))
  } # if(length(keep) < length(x.int))
  x.int = x.int[keep]
  y.int = y.int[keep]
  # Get initial sex estimates using mixture model clustering.
  mc = Mclust(cbind(x.int, y.int), G = 2, modelNames = "EEE")
  # The class with the higher X Chr means is female.
  f.clust = as.character(which.max(mc$parameters$mean[1,]))
  sex = mc$classification
  sex[sex == f.clust] = "F"
  sex[sex != "F"] = "M"
  # Use linear discriminant analysis to assign the sex.
  sex = factor(sex)
# DMG: lda stopped working with error "Error in t.default(La.res$vt) : 
# argument is not a matrix" on Nov. 28, 2013.
#  mod = lda(sex ~ x.int + y.int)
#  predict.sex = predict(mod)$class
# Using lm() instead.
  mod = lm(as.numeric(sex) ~ x.int + y.int)
  ps = predict(mod)
  predict.sex = rep("M", length(ps))
  predict.sex[ps < 1.5] = "F"
  predict.sex = factor(predict.sex)
  names(predict.sex) = names(sex)
  if(plot) {
    par(font = 2, font.axis = 2, font.lab = 2, las = 1)
    plot(x.int, y.int, pch = 16, col = as.numeric(predict.sex),
         xlab = "Mean X Chr Intensity", ylab = "Mean Y Chr Intensity",
         main = "Sex assignment")
    legend("topright", legend = levels(predict.sex), pch = 16, col = 1:2)
  } # if(plot)
  nm = names(predict.sex)
  predict.sex = as.character(predict.sex)
  names(predict.sex) = nm
  return(predict.sex)
} # sex.predict()
