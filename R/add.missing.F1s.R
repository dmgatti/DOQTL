################################################################################
# Add missing F1s to the X and Y matrices to act as pseudo-counts when we update
# genotype state means and variances.
# Arguments: founders: list containing founder information, either with 
#                      allele calls or intensities. Type is determined by the
#                      'method' attribute, which must be either 'allele' or 
#                      'intensity'.
#            snps: data.frame containing SNP information.
#            sampletype: character string indicating the type of 
#                        sample. One of c("DO", "CC", "DOF1", "HS", "other")
# Returns: list containing founder data, with the missing F1s added.
add.missing.F1s = function(founders, snps, sampletype = c("DO", "CC", "DOF1", 
                  "HS", "HSrat", "other")) {
  sampletype = match.arg(sampletype)
  missing = NULL
  if(any(nchar(founders$code) != 2)) {
    stop(paste("All founder codes must be 2 letters long. Please check",
	     "the founder codes and make sure that they are 2 letters long."))
  } # if(any(nchar(founders$code) != 2))
  
  if(length(unique(founders$sex)) == 2) {
    f.missing = founders$states$auto[which(!founders$states$auto %in% 
                founders$code[founders$sex == "F"])]
    f.missing.sex = rep("F", length(f.missing))
    m.missing = founders$states$auto[which(!founders$states$auto %in% 
                founders$code[founders$sex == "M"])]
    m.missing.sex = rep("M", length(m.missing))
    missing = data.frame(code = c(f.missing, m.missing), 
                         sex  = c(f.missing.sex, m.missing.sex))
    missing = missing[order(missing$code),]
  } else {
    if(unique(founders$sex) == "M") {
      if(attr(founders, "method") == "intensity") {
        warning("Only male founders found. X chromosome cannot be reconstructed.")
      } # if(attr(founders, "method") == "intensity")
      missing = founders$states$auto[which(!founders$states$auto %in% 
                founders$code[founders$sex == "M"])]
    } # if(unique(founders$sex) == "M")
    if(unique(founders$sex) == "F") {
      missing = founders$states$auto[which(!founders$states$auto %in% 
                founders$code[founders$sex == "F"])]
    } # if(unique(founders$sex) == "M")
  } # else
  mthd = attr(founders, "method")
  if(nrow(missing) > 0) {
    # Allele Calls
    if(mthd == "allele") {
      # Add columns for the new F1s.
      unique.geno = unique(unlist(apply(founders$geno, 2, unique)))
      stopifnot(all(unique.geno %in% c("A", "C", "G", "T", "H", "N")))
      founders$geno = data.frame(rbind(founders$geno, matrix("", nrow(missing), 
                  ncol(founders$geno), dimnames = list(missing$code, 
                  colnames(founders$geno)))),
                  row.names = make.unique(c(rownames(founders$geno), missing$code)))
      founders$code = c(founders$code, missing$code)
      names(founders$code) = rownames(founders$geno)
      founders$sex  = c(founders$sex, missing$sex)
      names(founders$sex) = rownames(founders$geno)
      founders$geno = lapply(founders$geno, factor, levels = unique.geno)
      founders$geno = data.frame(founders$geno)
      rownames(founders$geno) = names(founders$code)
      # Get the parents of each of the missing F1s.
      par = matrix(unlist(strsplit(missing$code, split = "")), nrow = 2)
      par = matrix(paste(par, par, sep = ""), nrow = 2)
      # For each F1, get the parental genotypes and figure out the F1 genotype.
      # Autosomes.
      curr.snps = which(snps[,2] %in% 1:19)
      # Get the parental genotypes.
      rn = sort(unique(as.vector(par)))
      par.gt = matrix("", length(rn), ncol(founders$geno), dimnames = 
               list(rn, colnames(founders$geno)))
      for(i in 1:nrow(par.gt)) {
        pg = lapply(founders$geno[founders$code == rownames(par.gt)[i],,
             drop = FALSE], table)
        par.gt[i,] = sapply(pg, function(z) { names(z)[which.max(z)] })
      } # for(i)
      # Fill in the F1 genotypes for the autosomes.
      for(i in 1:nrow(missing)) {
        geno.rows = which(founders$code == missing$code[i] & 
                          founders$sex  == missing$sex[i])
        p1g = par.gt[par[1,i],]
        p2g = par.gt[par[2,i],]
        homo = which(p1g == p2g)
        founders$geno[geno.rows, homo] = matrix(p1g[homo], nrow = length(geno.rows),
                                         ncol = length(homo), byrow = T)
        het = which(p1g != p2g)
        founders$geno[geno.rows, het] = "H"
      } # for(i)
      # X chromosome.
      curr.snps = which(snps[,2] == "X")
      if(length(curr.snps) > 0) {
        # Females
        missing.females = which(missing$sex == "F")
        for(i in missing.females) {
          geno.rows = which(founders$code == missing$code[i] & founders$sex == "F")
          p1g = par.gt[par[1,i],]
          p2g = par.gt[par[2,i],]
          homo = which(p1g == p2g)
          founders$geno[geno.rows, homo] = matrix(p1g[homo], nrow = length(geno.rows),
                                           ncol = length(homo), byrow = T)
          het = which(p1g != p2g)
          founders$geno[geno.rows, het] = "H"
        } # for(i)
        # Males
        # Because we don't know the maternal allele, set these to N. We don't
        # use the F1s to genotypes the males. We only use the founders.
        missing.males = which(missing$sex == "M")
        for(i in missing.males) {
          geno.rows = which(founders$code == missing$code[i] & founders$sex == "M")
          founders$geno[geno.rows, curr.snps] = "N"
        } # for(i)
      } # if(length(curr.snps) > 0)
      # Y chromosome.
      curr.snps = which(snps[,2] == "Y")
      if(length(curr.snps) > 0) {
        # Females
        # Set these to N.
        missing.females = which(missing$sex == "F")
        for(i in missing.females) {
          geno.rows = which(rownames(founders$geno) == missing$code[i] & 
                            founders$sex == "F")
          founders$geno[geno.rows, curr.snps] = "N"
        } # for(i)
        # Males
        # We don't used the F1s to genotype the Y chr, so set these to 'N'.
        missing.males = which(missing$sex == "M")
        for(i in missing.males) {
          geno.rows = which(rownames(founders$geno) == missing$code[i] & 
                            founders$sex == "F")
          founders$geno[geno.rows, curr.snps] = "N"
        } # for(i)
      } # if(length(curr.snps) > 0)
      # M chromosome.
      # We don't used the F1s to genotype the M chr, so set these to 'N'.
      curr.snps = which(snps[,2] == "M")
      if(length(curr.snps) > 0) {
        for(i in nrow(missing)) {
          geno.rows = which(rownames(founders$geno) == missing$code[i])
          founders$geno[geno.rows, curr.snps] = "N"
        } # for(i)
      } # if(length(curr.snps) > 0)    
      founders$geno = as.matrix(founders$geno)
    } else if(mthd == "intensity") {
      # Intensities
      # Add rows for the new F1s.
      missing$code = as.character(missing$code)
      founders$x = as.matrix(rbind(founders$x, matrix(0, length(missing$code), 
                   ncol(founders$x), dimnames = list(missing$code,
                   colnames(founders$x)))))
      founders$y = as.matrix(rbind(founders$y, matrix(0, length(missing$code), 
                   ncol(founders$y), dimnames = list(missing$code,
                   colnames(founders$y)))))
      founders$code = c(founders$code, missing$code)
      names(founders$code)[(length(founders$code) - 
                           nrow(missing) + 1):length(founders$code)] = missing$code
      founders$sex  = c(founders$sex, missing$sex)
      # Get the parents of each of the missing F1s.
      par = matrix(unlist(strsplit(missing$code, split = "")), nrow = 2)
      par = matrix(paste(par, par, sep = ""), nrow = 2)
      unique.par = unique(as.vector(par))
      par.means = list(x = matrix(0, length(unique.par), ncol(founders$x),
                       dimnames = list(unique.par, colnames(founders$x))), 
                       y = matrix(0, length(unique.par), ncol(founders$y),
                       dimnames = list(unique.par, colnames(founders$y))))
      for(i in 1:length(unique.par)) {
        rows = which(founders$code == unique.par[i])
        par.means$x[i,] = colMeans(founders$x[rows,,drop = FALSE], na.rm = TRUE)
        par.means$y[i,] = colMeans(founders$y[rows,,drop = FALSE], na.rm = TRUE)
      } # for(i)
      # Calculate the new F1s intensities as the mean of their parents.
      # Autosomes.
      curr.snps = which(snps[,2] %in% 1:19)
      for(i in seq(1, length(missing), 2)) {
        x.rows = which(rownames(founders$x) == missing$code[i])
        founders$x[x.rows, curr.snps] =
                   matrix(colMeans(par.means$x[par[,i], curr.snps], na.rm = TRUE),
                   length(x.rows), length(curr.snps), byrow = TRUE)
        founders$y[x.rows, curr.snps] = 
                   matrix(colMeans(par.means$y[par[,i], curr.snps], na.rm = TRUE),
                   length(x.rows), length(curr.snps), byrow = TRUE)
      } # for(i)
      # X chromosome.
      curr.snps = which(snps[,2] == "X")
      # Females
      missing.females = which(missing$sex == "F")
      for(i in missing.females) {
        x.rows = which(rownames(founders$x) == missing$code[i] & founders$sex == "F")
        founders$x[x.rows, curr.snps] = 
                   matrix(colMeans(par.means$x[par[,i], curr.snps], na.rm = TRUE),
                   length(x.rows), length(curr.snps), byrow = TRUE)
        founders$y[x.rows, curr.snps] = 
                   matrix(colMeans(par.means$y[par[,i], curr.snps], na.rm = TRUE),
                   length(x.rows), length(curr.snps), byrow = TRUE)
      } # for(i)
      # Males
      # We don't used the F1s for the male X chr, so set them to NA.
      missing.males = which(missing$sex == "M")
      for(i in missing.males) {
        x.rows = which(rownames(founders$x) == missing$code[i] & founders$sex == "M")
        founders$x[x.rows, curr.snps] = NA
        founders$y[x.rows, curr.snps] = NA
      } # for(i)
      # Y chromosome.
      curr.snps = which(snps[,2] == "Y")
      # Females
      # Females don't have a Y chr, so set them to NA.
      for(i in missing.females) {
        x.rows = which(rownames(founders$x) == missing$code[i] & founders$sex == "F")
        founders$x[x.rows, curr.snps] = NA
        founders$y[x.rows, curr.snps] = NA
      } # for(i)
      # Males
      # We don't used the F1s for the male Y chr, so set them to NA.      
      for(i in missing.males) {
        x.rows = which(rownames(founders$x) == missing$code[i] & founders$sex == "M")
        founders$x[x.rows, curr.snps] = NA
        founders$y[x.rows, curr.snps] = NA
      } # for(i)
      # M chromosome.
      curr.snps = which(snps[,2] == "M")
      # We don't used the F1s for the male Y chr, so set them to NA.      
      for(i in missing.males) {
        x.rows = which(rownames(founders$x) == missing$code[i] & founders$sex == "M")
        founders$x[x.rows, curr.snps] = NA
        founders$y[x.rows, curr.snps] = NA
      } # for(i)
      # Look for NaN in the founder data and try to replace them with
      # values.
      wh = which(colSums(is.nan(founders$x)) > 0 | colSums(is.nan(founders$y)) > 0)
      if(length(wh) > 0) {
        for(i in wh) {
          idx = which(is.nan(founders$x[,i]))
          for(j in idx) {
            cd = founders$code[j]
            stopifnot(length(unique(cd)) == 1)
            par = strsplit(cd, split = "")
            par = matrix(rep(unlist(par), each = 2), ncol = 2, byrow = TRUE)
            par = apply(par, 1, paste0, collapse = "")
            founders$x[j,i] = mean(founders$x[founders$code %in% par,i], na.rm = TRUE)
            founders$y[j,i] = mean(founders$y[founders$code %in% par,i], na.rm = TRUE)
          } # for(j)
        } # for(i)
      } # if(length(wh) > 0)
      if(any(is.nan(founders$x))) {
        stop("NaN in founders$x.")
      } # if(any(is.nan(founders$x)))
  
      if(any(is.nan(founders$y))) {
        stop("NaN in founders$y.")
      } # if(any(is.nan(founders$y)))
      if(sampletype == "DOF1") {
        keep = which(founders$code %in% founders$states$auto)
        founders$x = founders$x[keep,]
        founders$y = founders$y[keep,]
        founders$sex = founders$sex[keep]
        founders$code = founders$code[keep]
      } # if(sampletype == "DOF1")
    } else {
      stop(paste("Unknown method", mthd, "in add.missing.F1s()"))
    } # else
  } # if(length(missing) > 0)
  return(founders)
} # add.missing.F1s()
