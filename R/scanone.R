scanone = function(pheno, pheno.col = 1, probs = NULL, K = NULL, addcovar = NULL,
          intcovar = NULL, snps = NULL, model = c("additive", "full")) {

  model = match.arg(model)

  if(is.null(rownames(pheno))) {
    stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
  } # if(is.null(rownames(pheno)))

  if(is.null(addcovar)) {
    stop(paste("You must map using \\'sex\\' as an additive covariate. Also,",
         "we require sex to map on the X chromosome. We even require sex if",
         "you are only mapping with one sex."))
  } # if(missing(addcovar))

  sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
  if(length(sex.col) == 0) {
    stop(paste("addcovar must contain a column called \\'sex\\'. Please add",
         "a sex column with sex coded as either F & M or 0 for females and 1 for males."))
  } # if(length(sex.col) == 0)

  if(is.null(rownames(addcovar))) {
    stop("rownames(addcovar) is null. The sample IDs must be in rownames(addcovar).")
  } # if(is.null(rownames(addcovar)))

  if(!is.null(intcovar)) {
    if(is.null(rownames(intcovar))) {
      stop("rownames(intcovar) is null. The sample IDs must be in rownames(intcovar).")
    } # if(is.null(rownames(intcovar)))
  } # if(!is.null(intcovar))

  snps[,2] = as.character(snps[,2])

  num.auto = get.num.auto(snps)

  # Convert phenotype names to phenotype column numbers.
  if(is.character(pheno.col)) {
    pheno.col = match(pheno.col, colnames(pheno))
  } # if(is.character(pheno.col))

  # Get the intersection of all of the sample IDs from pheno, probs, K, 
  # addcovar and intcovar.
  tmp = synch.sample.IDs(pheno = pheno, probs = probs, K = K, addcovar = addcovar,
        intcovar = intcovar)
  pheno = tmp$pheno
  probs = tmp$probs
  K = tmp$K
  addcovar = tmp$addcovar
  intcovar = tmp$intcovar

  print(paste("Mapping with", nrow(pheno),"samples."))

  if(any(dim(probs) == 0)) {
    stop(paste("There are no matching samples in the data. Please",
         "verify that the sample IDs in rownames(pheno) match the sample",
         "IDs in rownames(probs), rownames(addcovar) and rownames(K)."))
  } # if(any(dim(probs) == 0))

  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]

  print(paste("Mapping with", nrow(snps), "markers."))

  if(any(dim(probs) == 0)) {
    stop(paste("There are no matching markers in snps and probs. Please",
         "verify that the marker IDs in snps[,1] match the marker",
         "IDs in dimnames(probs)[[3]]."))
  } # if(any(dim(probs) == 0))

  if(sum(rownames(pheno) %in% rownames(addcovar)) == 0) {
    stop(paste("rownames(pheno) does not contain any sample IDs in",
         "common with rownames(addcovar). Please make sure that the",
         "rownames in pheno and addcovar match."))
  } # if(sum(rownames(pheno) %in% rownames(addcovar)) == 0)

  addcovar = as.matrix(addcovar)

  # Match sample IDs in interactive covariates.
  if(!is.null(intcovar)) {

    intcovar = as.matrix(intcovar)
    intcovar = intcovar[rownames(intcovar) %in% rownames(pheno),,drop = FALSE]
    intcovar = intcovar[match(rownames(pheno), rownames(intcovar)),,drop = FALSE]
    if(is.null(colnames(intcovar))) {
      colnames(intcovar) = paste("intcovar", 1:ncol(intcovar), sep = ".")
    } # if(is.null(colnames(intcovar)))

  } # if(!is.null(intcovar))

  if(!is.null(K)) {

     # LOCO method.
     if(is.list(K)) {
       for(c in 1:length(K)) {
         K[[c]] = K[[c]][rownames(K[[c]]) %in% rownames(pheno), 
                         colnames(K[[c]]) %in% rownames(pheno)]
         K[[c]] = K[[c]][match(rownames(pheno), rownames(K[[c]])), 
                        match(rownames(pheno), colnames(K[[c]]))]
       } # for(c)
     } else {
       K = K[rownames(K) %in% rownames(pheno), colnames(K) %in% rownames(pheno)]
       K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), colnames(K))]
     } # else

  } # if(!is.null(K))

  retval = NULL
  if(is.null(K)) {
    retval = scanone.noK(pheno, pheno.col, probs, addcovar, intcovar, snps, model)
  } else if(is.list(K)) {
    retval = scanone.LOCO(pheno, pheno.col, probs, K, addcovar, intcovar, snps, model)
  } else {
    retval = scanone.K(pheno, pheno.col, probs, K, addcovar, intcovar, snps, model)
  } # else

  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)

  return(retval)

} # scanone()


################################################################################
# Help functions for scanone().
synch.sample.IDs = function(pheno = NULL, probs = NULL, K = NULL, addcovar = NULL,
                   intcovar = NULL) {

  samples = NULL

  # pheno
  if(!is.null(pheno)) {
      samples = rownames(pheno)
  } # if(!is.null(pheno))

  # probs
  if(!is.null(probs)) {

    if(is.null(samples)) {
      samples = rownames(probs)
    } else {
      samples = intersect(samples, rownames(probs))
    } # else

  } # if(!is.null(probs))

  # K
  if(!is.null(K)) {

    if(is.null(samples)) {

      if(is.list(K)) {
        samples = rownames(K[[1]])
      } else {
        samples = rownames(K)
      } # else

    } else {
      if(is.list(K)) {
        samples = intersect(samples, rownames(K[[1]]))
      } else {
        samples = intersect(samples, rownames(K))
      } # else
    } # else

  } # if(!is.null(probs))

  # addcovar
  if(!is.null(addcovar)) {

    if(is.null(samples)) {
      samples = rownames(addcovar)
    } else {
      samples = intersect(samples, rownames(addcovar))
    } # else

  } # if(!is.null(probs))

  # intcovar
  if(!is.null(intcovar)) {

    if(is.null(samples)) {
      samples = rownames(intcovar)
    } else {
      samples = intersect(samples, rownames(intcovar))
    } # else

  } # if(!is.null(probs))

  if(length(samples) == 0) {
    warning(paste("There were no samples in common among the variables",
            "passed in. Please make sure that there are rownames in",
            "common between all variables."))
  } # if(length(samples) == 0)

  # Now subset all of the variables and return them in a list.
  retval = NULL
  if(!is.null(pheno)) {
    pheno = pheno[samples,,drop = FALSE]
    retval = list(pheno = pheno)
  } # if(!is.null(pheno))

  if(!is.null(probs)) {
    probs = probs[samples,,,drop = FALSE]
    if(is.null(retval)) {
      retval = list(probs = probs)
    } else {
      retval[[length(retval) + 1]] = probs
      names(retval)[length(retval)] = "probs"
    } # else
  } # if(!is.null(probs))

  if(!is.null(K)) {

    if(is.list(K)) {
      for(i in 1:length(K)) {
        K[[i]] = K[[i]][samples, samples]
      } # for(i)
    } else {
      K = K[samples, samples]
    } # else

    if(is.null(retval)) {
      retval = list(K = K)
    } else {
      retval[[length(retval) + 1]] = K
      names(retval)[length(retval)] = "K"
    } # else
  } # if(!is.null(K))

  if(!is.null(addcovar)) {
    addcovar = addcovar[samples,,drop = FALSE]
    if(is.null(retval)) {
      retval = list(addcovar = addcovar)
    } else {
      retval[[length(retval) + 1]] = addcovar
      names(retval)[length(retval)] = "addcovar"
    } # else
  } # if(!is.null(addcovar))

  if(!is.null(intcovar)) {
    intcovar = intcovar[samples,,drop = FALSE]
    if(is.null(retval)) {
      retval = list(intcovar = intcovar)
    } else {
      retval[[length(retval) + 1]] = intcovar
      names(retval)[length(retval)] = "intcovar"
    } # else
  } # if(!is.null(intcovar))

  return(retval)

} # synch.sample.IDs()


scanone.noK = function(pheno, pheno.col, probs, addcovar, intcovar, snps, model) {

  num.auto = get.num.auto(snps)

  xchr = which(snps[,2] %in% "X")

  # We require sex to be in addcovar.
  if(length(xchr) > 0) {

    sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      stop(paste("You must map using \\'sex\\' as an additive covariate. Also,",
           "we require sex to map on the X chromosome. We even require sex if",
           "you are only mapping with one sex."))
    } # if(length(sex.col) == 0)
    addcovar[,sex.col] = as.numeric(factor(addcovar[,sex.col])) - 1

  } # if(length(xchr) > 0)

  # Make the results list.
  retval = as.list(1:length(pheno.col))
  names(retval) = colnames(pheno)[pheno.col]
  index = 1

  for(i in pheno.col) {

    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p) & !is.nan(p) & !is.infinite(p))

    # Autosomes
    auto.qtl = NULL
    if(!is.na(num.auto)) {

      auto = which(snps[,2] %in% 1:num.auto)

      # With covariates.
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                   rowSums(is.nan(addcovar)) == 0 &
                                   rowSums(is.infinite(addcovar)) == 0))

      if(is.null(intcovar)) {

        # Additive covariates only.
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                   addcovar = addcovar[keep,,drop = FALSE],
                   snps = snps[auto,])

      } else {

        auto.qtl = qtl.qtlrel(pheno = p[keep], probs = probs[keep,,auto],
                   addcovar = addcovar[keep,,drop = FALSE], 
                   intcovar = intcovar[keep,,drop = FALSE], 
                   snps = snps[auto,])

      } # else

      auto.qtl = list(lod = list(A = auto.qtl$lod), 
                      coef = list(A = auto.qtl$coef))

    } # if(!is.na(num.auto))

    # X chromosome.
    if(length(xchr) > 0) {

      # Get the sex from addcovar. We forced it to be 0 or 1 above.
      sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      females = which(addcovar[,sex.col] == 0)
      males   = which(addcovar[,sex.col] == 1)
      mfprobs = NULL

      if(length(females) > 0 & length(males) > 0) {

        # Take the male and female probabilities and place them into 
        # one big array.
        if(model == "additive") {

          mfprobs = array(0, c(dim(probs)[1], 2 * dim(probs)[2], length(xchr)),
                    dimnames = list(dimnames(probs)[[1]], paste(rep(c("F", "M"), 
                    each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."), 
                    dimnames(probs)[[3]][xchr]))

          for(j in females) {

            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]

          } # for(j)

          for(j in males) {

            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,,xchr]

          } # for(j)

        } else if(model == "full") {

          tmp = matrix(unlist(strsplit(dimnames(probs)[[2]], split = "")),
                nrow = 2)
          homo = dimnames(probs)[[2]][which(tmp[1,] == tmp[2,])]
          mfprobs = array(0, c(dim(probs)[1], dim(probs)[2] + length(homo),
                    length(xchr)),
                    dimnames = list(dimnames(probs)[[1]], c(paste(rep("F", 
                    each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."),
                    paste("M", homo, sep = ".")), dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,homo,xchr]
          } # for(j)

        } # else if(model == "full")

        # If we have both males and females, then we need to remove one 
        # column from the males.
        mfprobs = mfprobs[,-grep("M.A", dimnames(mfprobs)[[2]]),]

      } else if(length(females) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("F", dimnames(mfprobs)[[2]], sep = ".")
      } else if(length(males) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("M", dimnames(mfprobs)[[2]], sep = ".")
      } # else if(length(males) > 0)


      # With covariates.
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))

      if(is.null(intcovar)) {

        # Additive covariates only.
        x.qtl = fast.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                addcovar = addcovar[keep,-1,drop = FALSE],
                snps = snps[xchr,])

      } else {
        # Additive & interactive covariates.
        x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                addcovar = addcovar[keep,,drop = FALSE], 
                intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])

      } # else

      if(!is.null(auto.qtl)) {

        auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
        auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)

      } else {

        auto.qtl = x.qtl

      } # else

    } # if(length(xchr) > 0)

    retval[[index]] = auto.qtl
    class(retval[[index]]) = c("doqtl", class(retval[[index]]))
    attr(retval[[index]], "model") = "additive"
    index = index + 1

  } # for(i)

  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)

  return(retval)

} # scanone.noK()


################################################################################
scanone.K = function(pheno, pheno.col = 1, probs, K, addcovar, intcovar, snps, model) {

  num.auto = get.num.auto(snps)

  xchr = which(snps[,2] %in% "X")

  # We require sex to be in addcovar.
  sex = NULL
  if(length(xchr) > 0) {

    sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      stop(paste("You must map using \\'sex\\' as an additive covariate. Also,",
           "we require sex to map on the X chromosome. We even require sex if",
           "you are only mapping with one sex."))
    } # if(length(sex.col) == 0)
    addcovar[,sex.col] = as.numeric(factor(addcovar[,sex.col])) - 1

  } # if(length(xchr) > 0)

  # Make the results list.
  retval = as.list(1:length(pheno.col))
  names(retval) = colnames(pheno)[pheno.col]
  index = 1
  for(i in pheno.col) {

    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p) & !is.nan(p) & !is.infinite(p))

    # Autosomes
    auto.qtl = NULL
    if(!is.na(num.auto)) {

      auto = which(snps[,2] %in% 1:num.auto)

      # With covariates.
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                   rowSums(is.nan(addcovar)) == 0 &
                                   rowSums(is.infinite(addcovar)) == 0))
      if(is.null(intcovar)) {

        # Additive covariates only.
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                   K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                   snps = snps[auto,])

      } else {

        auto.qtl = qtl.qtlrel(pheno = p[keep], probs = probs[keep,,auto],
                   K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                   intcovar = intcovar[keep,,drop = FALSE], snps = snps[auto,])

      } # else

      auto.qtl = list(lod = list(A = auto.qtl$lod), 
                      coef = list(A = auto.qtl$coef))

    } # if(!is.na(num.auto))

    # X chromosome.
    if(length(xchr) > 0) {

      # Get the sex from addcovar. We forced it to be 0 or 1 above.
      sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      females = which(addcovar[,sex.col] == 0)
      males   = which(addcovar[,sex.col] == 1)
      mfprobs = NULL

      if(length(females) > 0 & length(males) > 0) {

        # Take the male and female probabilities and place them into 
        # one big array.
        if(model == "additive") {

          mfprobs = array(0, c(dim(probs)[1], 2 * dim(probs)[2], length(xchr)),
                    dimnames = list(dimnames(probs)[[1]], paste(rep(c("F", "M"), 
                    each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."), 
                    dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,,xchr]
          } # for(j)

        } else if(model == "full") {

          tmp = matrix(unlist(strsplit(dimnames(probs)[[2]], split = "")),
                nrow = 2)
          homo = dimnames(probs)[[2]][which(tmp[1,] == tmp[2,])]
          mfprobs = array(0, c(dim(probs)[1], dim(probs)[2] + length(homo),
                    length(xchr)),
                    dimnames = list(dimnames(probs)[[1]], c(paste(rep("F", 
                    each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."),
                    paste("M", homo, sep = ".")), dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,homo,xchr]
          } # for(j)

        } # else if(model == "full")

        # If we have both males and females, then we need to remove one 
        # column from the males.
        mfprobs = mfprobs[,-grep("M.A", dimnames(mfprobs)[[2]]),]

      } else if(length(females) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("F", dimnames(mfprobs)[[2]], sep = ".")
      } else if(length(males) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("M", dimnames(mfprobs)[[2]], sep = ".")
      } # else if(length(males) > 0)

      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
      if(is.null(intcovar)) {

        # Additive covariates only.
        x.qtl = fast.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                snps = snps[xchr,])

      } else {

        # Additive & interactive covariates.
        x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])

      } # else

      if(!is.null(auto.qtl)) {
        auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
        auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)
      } else {
        auto.qtl = x.qtl
      } # else

    } # if(length(xchr) > 0)

    retval[[index]] = auto.qtl
    class(retval[[index]]) = c("doqtl", class(retval[[index]]))
    attr(retval[[index]], "model") = "additive"
    index = index + 1

  } # for(i)

  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)

  return(retval)

} # scanone.K()


################################################################################
scanone.LOCO = function(pheno, pheno.col = 1, probs, K, addcovar, intcovar, snps, model) {

  num.auto = get.num.auto(snps)

  xchr = which(snps[,2] %in% "X")

  # We require sex to be in addcovar.
  sex = NULL
  if(length(xchr) > 0) {

    sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      stop(paste("You must map using \\'sex\\' as an additive covariate. Also,",
           "we require sex to map on the X chromosome. We even require sex if",
           "you are only mapping with one sex."))
    } # if(length(sex.col) == 0)
    addcovar[,sex.col] = as.numeric(factor(addcovar[,sex.col])) - 1

  } # if(length(xchr) > 0)

  # Make the results list.
  retval = as.list(1:length(pheno.col))
  names(retval) = colnames(pheno)[pheno.col]
  index = 1
  for(i in pheno.col) {

    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p) & !is.nan(p) & !is.infinite(p))

    # Autosomes
    auto.qtl = NULL
    if(!is.na(num.auto)) {

      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                   rowSums(is.nan(addcovar)) == 0 &
                                   rowSums(is.infinite(addcovar)) == 0))

      if(is.null(intcovar)) {

        # Additive covariates only.
        snprng = which(snps[,2] == 1)
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,snprng], 
                   K = K[[1]][keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                   snps = snps[snprng,])

        for(c in 2:num.auto) {
          snprng = which(snps[,2] == c)
          tmp = fast.qtlrel(pheno = p[keep, drop = FALSE], probs = probs[keep,,snprng], 
                K = K[[c]][keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                snps = snps[snprng,])
          auto.qtl$lod  = rbind(auto.qtl$lod,  tmp$lod)
          auto.qtl$coef = rbind(auto.qtl$coef, tmp$coef)
        } # for(c)

      } else {

        # Additive and interactive covariates.
        snprng = which(snps[,2] == 1)
        auto.qtl = qtl.qtlrel(pheno = p[keep, drop = FALSE], 
                   probs = probs[keep,,snprng], K = K[[1]][keep,keep], 
                   addcovar = addcovar[keep,,drop = FALSE], 
                   intcovar = intcovar[keep,,drop = FALSE], snps = snps[snprng,])
        for(c in 2:num.auto) {
          snprng = which(snps[,2] == c)
          tmp = qtl.qtlrel(pheno = p[keep, drop = FALSE],
                probs = probs[keep,,snprng], K = K[[c]][keep,keep], 
                addcovar = addcovar[keep,,drop = FALSE], 
                intcovar = intcovar[keep,,drop = FALSE], snps = snps[snprng,])
           auto.qtl$lod  = rbind(auto.qtl$lod,  tmp$lod)
           auto.qtl$coef = rbind(auto.qtl$coef, tmp$coef)
        } # for(c)

      } # else

      auto.qtl = list(lod  = list(A = auto.qtl$lod), 
                      coef = list(A = auto.qtl$coef))

    } # if(!is.na(num.auto))

    # X chromosome.
    if(length(xchr) > 0) {

      # Get the sex from addcovar. We forced it to be 0 or 1 above.
      sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      females = which(addcovar[,sex.col] == 0)
      males   = which(addcovar[,sex.col] == 1)
      mfprobs = NULL

      if(length(females) > 0 & length(males) > 0) {

        # Take the male and female probabilities and place them into 
        # one big array.
        if(model == "additive") {

          mfprobs = array(0, c(nrow(probs), 2 * ncol(probs), length(xchr)),
                    dimnames = list(dimnames(probs)[[1]], paste(rep(c("F", "M"), 
                    each = ncol(probs)), dimnames(probs)[[2]], sep = "."), 
                    dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,,xchr]
          } # for(j)

        } else if(model == "full") {

          tmp = matrix(unlist(strsplit(dimnames(probs)[[2]], split = "")),
                nrow = 2)
          homo = dimnames(probs)[[2]][which(tmp[1,] == tmp[2,])]
          mfprobs = array(0, c(dim(probs)[1], dim(probs)[2] + length(homo),
                    length(xchr)),
                    dimnames = list(dimnames(probs)[[1]], c(paste(rep("F", 
                    each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."),
                    paste("M", homo, sep = ".")), dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,homo,xchr]
          } # for(j)

        } # else if(model == "full")

        # If we have both males and females, then we need to remove one 
        # column from the males.
        mfprobs = mfprobs[,-grep("M.A", dimnames(mfprobs)[[2]]),]

      } else if(length(females) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("F", dimnames(mfprobs)[[2]], sep = ".")
      } else if(length(males) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("M", dimnames(mfprobs)[[2]], sep = ".")
      } # else if(length(males) > 0)

      x.qtl = NULL
      if(is.null(intcovar)) {

        # Additive covariates only.
        x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                K = K[["X"]][keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                snps = snps[xchr,])

      } else {

        # Additive & interactive covariates.
        x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                K = K[["X"]][keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])

      } # else

      if(!is.null(auto.qtl)) {
        auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
        auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)
      } else {
        auto.qtl = x.qtl
      } # else

    } # if(length(xchr) > 0)

    retval[[index]] = auto.qtl
    class(retval[[index]]) = c("doqtl", class(retval[[index]]))
    attr(retval[[index]], "model") = "additive"
    index = index + 1

  } # for(i)

  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)

  return(retval)

} # scanone.LOCO()

