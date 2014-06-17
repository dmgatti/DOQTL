################################################################################
# Add missing F1s to the X and Y matrices to act as pseudo-counts when we update
# genotype state means and variances.
# Arguments: founders: list containing founder information, either with 
#                      allele calls or intensities. Type is determined by the
#                      'method' attribute, which must be either 'allele' or 
#                      'intensity'.
#            snps: data.frame containing SNP information.
# Returns: list containing founder data, with the missing F1s added.
add.missing.F1s = function(founders, snps) {
  missing = founders$states$auto[which(!founders$states$auto %in% 
            founders$code[founders$sex == "F"])]
  mthd = attr(founders, "method")
  if(length(missing) > 0) {
    # Expand missing to 2 samples per genotype, one for males and one for females.
    missing = rep(missing, each = 2)
    # Allele Calls
    if(mthd == "allele") {
      # Add columns for the new F1s.
      founders$geno = data.frame(rbind(founders$geno, matrix("", length(missing), 
                      ncol(founders$geno),dimnames = list(missing, 
                      colnames(founders$geno)))), check.names = FALSE)
      founders$code = c(founders$code, missing)
      names(founders$code)[(length(founders$code) - 
                           length(missing) + 1):length(founders$code)] = missing
      founders$sex  = c(founders$sex, rep(c("F", "M"), length(missing) / 2))
      names(founders$sex)[(length(founders$sex) - 
                           length(missing) + 1):length(founders$sex)] = missing
      unique.geno = unique(unlist(apply(founders$geno, 2, unique)))
      founders$geno = lapply(founders$geno, factor, levels = unique.geno)
      founders$geno = data.frame(founders$geno)
      # Get the parents of each of the missing F1s.
      par = matrix(unlist(strsplit(missing, split = "")), nrow = 2)
      par = matrix(paste(par, par, sep = ""), nrow = 2)
      # For each F1, get the parental genotypes and figure out the F1 genotype.
      # Autosomes.
      curr.snps = which(snps[,2] %in% 1:19)
      for(i in seq(1, length(missing), 2)) {
        geno.rows = which(founders$code == missing[i])
        p1g = lapply(founders$geno[founders$code == par[1,i], curr.snps, drop = FALSE],
              table)
        p2g = lapply(founders$geno[founders$code == par[2,i], curr.snps, drop = FALSE],
              table)
        p1g = sapply(p1g, function(z) { names(z)[which.max(z)] })
        p2g = sapply(p2g, function(z) { names(z)[which.max(z)] })
        homo = which(p1g == p2g)
        founders$geno[geno.rows, homo] = p1g[homo]
        het = which(p1g != p2g)
        founders$geno[geno.rows, het] = "H"
      } # for(i)
      # X chromosome.
      curr.snps = which(snps[,2] == "X")
      if(length(curr.snps) > 0) {
        # Females
        for(i in seq(1, length(missing), 2)) {
          geno.rows = which(founders$code == missing[i] & founders$sex == "F")
          p1g = lapply(founders$geno[founders$code == par[1,i], curr.snps, drop = FALSE],
                table)
          p2g = lapply(founders$geno[founders$code == par[2,i], curr.snps, drop = FALSE],
                table)
          p1g = sapply(p1g, function(z) { names(z)[which.max(z)] })
          p2g = sapply(p2g, function(z) { names(z)[which.max(z)] })
          homo = which(p1g == p2g)
          founders$geno[geno.rows, curr.snps][,homo] = p1g[homo]
          het = which(p1g != p2g)
          founders$geno[geno.rows, curr.snps][,het] = "H"
        } # for(i)
            
        # Males
        # Because we don't know the maternal allele, set these to N. We don't
        # use the F1s to genotypes the males. We only use the founders.
        for(i in seq(1, length(missing), 2)) {
          geno.rows = which(founders$code == missing[i] & founders$sex == "M")
          founders$geno[geno.rows, curr.snps, drop = FALSE] = "N"
        } # for(i)
      } # if(length(curr.snps) > 0)
      # Y chromosome.
      curr.snps = which(snps[,2] == "Y")
      if(length(curr.snps) > 0) {
        # Females
        # Set these to N.
        for(i in seq(1, length(missing), 2)) {
          geno.rows = which(rownames(founders$geno) == missing[i] & founders$sex == "F")
          founders$geno[geno.rows, curr.snps, drop = FALSE] = "N"
        } # for(i)
        # Males
        # We don't used the F1s to genotype the Y chr, so set these to 'N'.
        for(i in seq(1, length(missing), 2)) {
          geno.rows = which(rownames(founders$geno) == missing[i] & founders$sex == "F")
          founders$geno[geno.rows, curr.snps, drop = FALSE] = "N"
        } # for(i)
      } # if(length(curr.snps) > 0)
      # M chromosome.
      # We don't used the F1s to genotype the M chr, so set these to 'N'.
      curr.snps = which(snps[,2] == "M")
      if(length(curr.snps) > 0) {
        for(i in seq(1, length(missing), 2)) {
          geno.rows = which(rownames(founders$geno) == missing[i] & founders$sex == "F")
          founders$geno[geno.rows, curr.snps, drop = FALSE] = "N"
        } # for(i)
      } # if(length(curr.snps) > 0)    
      founders$geno = as.matrix(founders$geno)
    } else if(mthd == "intensity") {
    # Intensities
  
      # Add rows for the new F1s.
      founders$x = as.matrix(rbind(founders$x, matrix(0, length(missing), 
                   ncol(founders$x), dimnames = list(missing,
                   colnames(founders$x)))))
      founders$y = as.matrix(rbind(founders$y, matrix(0, length(missing), 
                   ncol(founders$y), dimnames = list(missing,
                   colnames(founders$y)))))
      founders$code = c(founders$code, missing)
      names(founders$code)[(length(founders$code) - 
                           length(missing) + 1):length(founders$code)] = missing
      founders$sex  = c(founders$sex, rep(c("F", "M"), length(missing) / 2))
      names(founders$sex)[(length(founders$sex) - 
                           length(missing) + 1):length(founders$sex)] = missing
      # Get the parents of each of the missing F1s.
      par = matrix(unlist(strsplit(missing, split = "")), nrow = 2)
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
        x.rows = which(rownames(founders$x) == missing[i])
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
      for(i in seq(1, length(missing), 2)) {
        x.rows = which(rownames(founders$x) == missing[i] & founders$sex == "F")
        founders$x[x.rows, curr.snps] = 
                   matrix(colMeans(par.means$x[par[,i], curr.snps], na.rm = TRUE),
                   length(x.rows), length(curr.snps), byrow = TRUE)
        founders$y[x.rows, curr.snps] = 
                   matrix(colMeans(par.means$y[par[,i], curr.snps], na.rm = TRUE),
                   length(x.rows), length(curr.snps), byrow = TRUE)
      } # for(i)
      # Males
      # We don't used the F1s for the male X chr, so set them to NA.
      for(i in seq(1, length(missing), 2)) {
        x.rows = which(rownames(founders$x) == missing[i] & founders$sex == "M")
        founders$x[x.rows, curr.snps] = NA
        founders$y[x.rows, curr.snps] = NA
      } # for(i)
      # Y chromosome.
      curr.snps = which(snps[,2] == "Y")
      # Females
      # Females don't have a Y chr, so set them to NA.
      for(i in seq(1, length(missing), 2)) {
        x.rows = which(rownames(founders$x) == missing[i] & founders$sex == "F")
        founders$x[x.rows, curr.snps] = NA
        founders$y[x.rows, curr.snps] = NA
      } # for(i)
      # Males
      # We don't used the F1s for the male Y chr, so set them to NA.      
      for(i in seq(1, length(missing), 2)) {
        x.rows = which(rownames(founders$x) == missing[i] & founders$sex == "M")
        founders$x[x.rows, curr.snps] = NA
        founders$y[x.rows, curr.snps] = NA
      } # for(i)
      # M chromosome.
      curr.snps = which(snps[,2] == "M")
      # We don't used the F1s for the male Y chr, so set them to NA.      
      for(i in seq(1, length(missing), 2)) {
        x.rows = which(rownames(founders$x) == missing[i] & founders$sex == "M")
        founders$x[x.rows, curr.snps] = NA
        founders$y[x.rows, curr.snps] = NA
      } # for(i)
    } else {
      stop(paste("Unknown method", mthd, "in add.missing.F1s()"))
    } # else
  } # if(length(missing) > 0)
  return(founders)
} # add.missing.F1s()
