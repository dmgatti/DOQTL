scanone = function(pheno, pheno.col = 1, probs, K, addcovar, intcovar, snps,
          model = c("additive", "full")) {

  model = match.arg(model)

  if(is.null(rownames(pheno))) {
    stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
  } # if(is.null(rownames(pheno)))

  if(!missing(addcovar)) {
    if(is.null(rownames(addcovar))) {
      stop("rownames(addcovar) is null. The sample IDs must be in rownames(addcovar).")
    } # if(is.null(rownames(addcovar)))
  } # if(!missing(addcovar))

  if(!missing(intcovar)) {
    if(is.null(rownames(intcovar))) {
      stop("rownames(intcovar) is null. The sample IDs must be in rownames(intcovar).")
    } # if(is.null(rownames(intcovar)))
  } # if(!missing(intcovar))

  snps[,2] = as.character(snps[,2])

  num.auto = get.num.auto(snps)

  # Convert phenotype names to phenotype column numbers.
  if(is.character(pheno.col)) {
    pheno.col = match(pheno.col, colnames(pheno))
  } # if(is.character(pheno.col))

  # Get the intersection of all of the sample IDs from pheno, probs, K, 
  # addcovar and intcovar.
  samples = intersect(rownames(pheno), rownames(probs))
  if(!missing(K)) {
    if(is.list(K)) {
      # We're assuming that all K matrices have the same rownames.
      samples = intersect(samples, rownames(K[[1]]))
    } else {
      samples = intersect(samples, rownames(K))
    } # else
  } # if(!missing(K))
  if(!missing(addcovar)) {
    samples = intersect(samples, rownames(addcovar))
  } # if(!missing(addcovar))
  if(!missing(intcovar)) {
    samples = intersect(samples, rownames(intcovar))
  } # if(!missing(intcovar))

  # Subset the data to include only samples in common.
  pheno = pheno[samples,, drop = FALSE]
  probs = probs[samples,,]

  if(!missing(K)) {
    if(is.list(K)) {
      for(c in 1:length(K)) {
        K[[c]] = K[[c]][samples, samples]
      } # for(c)
    } else {
        K = K[samples, samples]
    } # else
  } # if(!missing(K))

  if(!missing(addcovar)) {
    addcovar = addcovar[samples,, drop = FALSE]
  } # if(!missing(addcovar))

  if(!missing(intcovar)) {
    intcovar = intcovar[samples,, drop = FALSE]
  } # if(!missing(intcovar))

  message(paste("Mapping with", nrow(pheno),"samples.\n"))

  if(any(dim(probs) == 0)) {
    stop(paste("There are no matching samples in pheno and probs. Please",
         "verify that the sample IDs in rownames(pheno) match the sample",
         "IDs in dimnames(probs)[[1]]."))
  } # if(any(dim(probs) == 0))

#  probs = filter.geno.probs(probs)
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]

  # Match sample IDs in additive covariates.
  if(!missing(addcovar)) {

    addcovar = data.frame(addcovar, stringsAsFactors = TRUE)
    rn = rownames(addcovar)
    addcovar = lapply(addcovar, as.numeric)
    cn = names(addcovar)
    addcovar = matrix(unlist(addcovar), length(addcovar[[1]]),
               length(addcovar), dimnames = list(rn, cn))
    rownames(addcovar) = rn
    addcovar = addcovar[rownames(addcovar) %in% rownames(pheno),,drop = FALSE]
    addcovar = addcovar[match(rownames(pheno), rownames(addcovar)),,drop = FALSE]

    if(sum(rownames(pheno) %in% rownames(addcovar)) == 0) {
      stop(paste("rownames(pheno) does not contain any sample IDs in",
           "common with rownames(addcovar). Please make sure that the",
           "rownames in pheno and addcovar match."))
    } # if(sum(rownames(pheno) %in% rownames(addcovar)) == 0)

    pheno = pheno[rownames(pheno) %in% rownames(addcovar),,drop = FALSE]
    addcovar = addcovar[rownames(addcovar) %in% rownames(pheno),,drop = FALSE]
    addcovar = addcovar[match(rownames(pheno), rownames(addcovar)),,drop = FALSE]
    probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]

    if(is.null(colnames(addcovar))) {
      colnames(addcovar) = paste("addcovar", 1:ncol(addcovar), sep = ".")
    } # if(is.null(colnames(addcovar)))

    addcovar = as.matrix(addcovar)
  } # if(!is.null(addcovar))

  # Match sample IDs in interactive covariates.
  if(!missing(intcovar)) {

    intcovar = as.matrix(intcovar)
    intcovar = intcovar[rownames(intcovar) %in% rownames(pheno),,drop = FALSE]
    intcovar = intcovar[match(rownames(pheno), rownames(intcovar)),,drop = FALSE]
    if(is.null(colnames(intcovar))) {
      colnames(intcovar) = paste("intcovar", 1:ncol(intcovar), sep = ".")
    } # if(is.null(colnames(intcovar)))

  } # if(!is.null(intcovar))

  if(!missing(K)) {

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

  } # if(!missing(K))

  retval = NULL
  if(missing(K)) {
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
synch.sample.IDs = function(vars) {

  use = which(sapply(vars, is.null))
  samples = unique(unlist(lapply(vars[use], rownames)))
  for(i in use) {
    if(length(dim(vars[[i]])) == 2) {
      vars[[i]] = vars[[i]][samples,,drop = FALSE]
    } else if(length(dim(vars[[i]])) == 3) {
      vars[[i]] = vars[[i]][samples,,,drop = FALSE]
    }
  } # for(i)

  return(vars)

} # synch.sample.IDs()


scanone.noK = function(pheno, pheno.col, probs, addcovar, intcovar, snps, model) {

  num.auto = get.num.auto(snps)

  xchr = which(snps[,2] %in% "X")
  sex = NULL
  if(length(xchr) > 0) {
    sex.col = grep("^sex$", colnames(pheno), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      if(!missing(addcovar)) {
        sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      } # if(!missing(addcovar))

      if(length(sex.col) == 0) {
        stop(paste("There is no sex column in the phenotypes or covariates.",
                   "Please add a sex column for proper mapping on the X chromosome."))
      } else {
        sex = addcovar[,sex.col]
      } # else

    } else {
      sex = pheno[,sex.col]
    } # else
    sex = toupper(sex)
    # This takes care of the case where "F" gets turned into "FALSE",
    # by R.
    wh = which(sex == "FALSE")
    if(length(wh) > 0) {
      sex[wh] = "F"
    } # if(is.logical(sex))
    names(sex) = rownames(pheno)

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
    auto = which(snps[,2] %in% 1:num.auto)

    if(missing(addcovar)) {

      # No covariates.
      auto.qtl = fast.qtlrel(pheno = p[keep,drop = FALSE], probs = probs[keep,,auto], 
                 snps = snps[auto,])

    } else {

      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                   rowSums(is.nan(addcovar)) == 0 &
                                   rowSums(is.infinite(addcovar)) == 0))
      if(missing(intcovar)) {

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
    } # else

    auto.qtl = list(lod = list(A = auto.qtl$lod), 
                    coef = list(A = auto.qtl$coef))

    # X chromosome.
    xchr = which(snps[,2] %in% "X")

    if(length(xchr) > 0) {

      females = which(sex == "F" | sex == "0")
      males   = which(sex == "M" | sex == "1")
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

      if(missing(addcovar)) {

        # No covariates, but we need sex for the null hypothesis.
        if(length(grep("^sex$", colnames(pheno), ignore.case = TRUE)) > 0) {
          sex = pheno[,grep("^sex$", colnames(pheno), ignore.case = TRUE)]
          sex = as.matrix(as.numeric(sex == "M"))
          dimnames(sex) = list(rownames(pheno), "sex")
        } # if(length(grep("^sex$", colnames(pheno), ignore.case = TRUE)) > 0)

        x.qtl = qtl.qtlrel(pheno = p[keep,drop = FALSE], probs = mfprobs[keep,,],
                snps = snps[xchr,,drop = FALSE])

      } else {

        keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
        if(missing(intcovar)) {

          # Additive covariates only.
          x.qtl = fast.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                  addcovar = addcovar[keep,,drop = FALSE],
                  snps = snps[xchr,])

        } else {
          # Additive & interactive covariates.
          x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                  addcovar = addcovar[keep,,drop = FALSE], 
                  intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])

        } # else
      } # else

      auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
      auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)

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
  sex = NULL
  if(length(xchr) > 0) {
    sex.col = grep("^sex$", colnames(pheno), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      if(!missing(addcovar)) {
        sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      } # if(!missing(addcovar))

      if(length(sex.col) == 0) {
        stop(paste("There is no sex column in the phenotypes or covariates.",
                   "Please add a sex column for proper mapping on the X chromosome."))
      } else {
        sex = addcovar[,sex.col]
      } # else

    } else {
      sex = pheno[,sex.col]
    } # else
    sex = toupper(sex)
    # This takes care of the case where "F" gets turned into "FALSE",
    # by R.
    wh = which(sex == "FALSE")
    if(length(wh) > 0) {
      sex[wh] = "F"
    } # if(is.logical(sex))
    names(sex) = rownames(pheno)

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
    auto = which(snps[,2] %in% 1:num.auto)

    if(missing(addcovar)) {

      # No covariates.
      auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                 K = K[keep,keep], snps = snps[auto,])

    } else {

      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                   rowSums(is.nan(addcovar)) == 0 &
                                   rowSums(is.infinite(addcovar)) == 0))
      if(missing(intcovar)) {

        # Additive covariates only.
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                   K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                   snps = snps[auto,])

      } else {

        auto.qtl = qtl.qtlrel(pheno = p[keep], probs = probs[keep,,auto],
                   K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                   intcovar = intcovar[keep,,drop = FALSE], snps = snps[auto,])

      } # else
    } # else

    auto.qtl = list(lod = list(A = auto.qtl$lod), 
                    coef = list(A = auto.qtl$coef))

    # X chromosome.
    xchr = which(snps[,2] %in% "X")

    if(length(xchr) > 0) {

      females = which(sex == "F" | sex == "0")
      males   = which(sex == "M" | sex == "1")
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

      if(missing(addcovar)) {

        # No covariates, but we need sex for the null hypothesis.
        if(length(grep("^sex$", colnames(pheno), ignore.case = TRUE)) > 0) {
          sex = pheno[,grep("^sex$", colnames(pheno), ignore.case = TRUE)]
          sex = as.matrix(as.numeric(sex == "M"))
          dimnames(sex) = list(rownames(pheno), "sex")
        } # if(length(grep("^sex$", colnames(pheno), ignore.case = TRUE)) > 0)

        x.qtl = qtl.qtlrel(pheno = p[keep,drop = FALSE], probs = mfprobs[keep,,],
                K = K[keep,keep], snps = snps[xchr,,drop = FALSE])

      } else {

        keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
        if(missing(intcovar)) {

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
      } # else

      auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
      auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)

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
  sex = NULL
  if(length(xchr) > 0) {
    sex.col = grep("^sex$", colnames(pheno), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      if(!missing(addcovar)) {
        sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      } # if(!missing(addcovar))

      if(length(sex.col) == 0) {
        stop(paste("There is no sex column in the phenotypes or covariates.",
                   "Please add a sex column for proper mapping on the X chromosome."))
      } else {
        sex = addcovar[,sex.col]
      } # else

    } else {
      sex = pheno[,sex.col]
    } # else
    sex = toupper(sex)
    # This takes care of the case where "F" gets turned into "FALSE",
    # by R.
    wh = which(sex == "FALSE")
    if(length(wh) > 0) {
      sex[wh] = "F"
    } # if(is.logical(sex))
    names(sex) = rownames(pheno)

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
    if(missing(addcovar)) {

      # No covariates.
      snprng = which(snps[,2] == 1)
      auto.qtl = fast.qtlrel(pheno = p[keep, drop = FALSE], 
                 probs = probs[keep,,snprng], 
                 K = K[[1]][keep,keep], snps = snps[snprng,])
      for(c in 2:num.auto) {
        snprng = which(snps[,2] == c)
        tmp = fast.qtlrel(pheno = p[keep, drop = FALSE],
              probs = probs[keep,,snprng], 
              K = K[[c]][keep,keep], snps = snps[snprng,])
        auto.qtl$lod  = rbind(auto.qtl$lod, tmp$lod)
        auto.qtl$coef = rbind(auto.qtl$coef, tmp$coef)
      } # for(c)

    } else {

      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                   rowSums(is.nan(addcovar)) == 0 &
                                   rowSums(is.infinite(addcovar)) == 0))

      if(missing(intcovar)) {

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
    } # else

    auto.qtl = list(lod  = list(A = auto.qtl$lod), 
                    coef = list(A = auto.qtl$coef))

    # X chromosome.
    xchr = which(snps[,2] %in% "X")

    if(length(xchr) > 0) {

      females = which(sex == "F" | sex == "0")
      males   = which(sex == "M" | sex == "1")
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

      if(missing(addcovar)) {

        # No covariates, but we need sex for the null hypothesis.
        if(length(grep("^sex$", colnames(pheno), ignore.case = TRUE)) > 0) {
          sex = pheno[,grep("^sex$", colnames(pheno), ignore.case = TRUE)]
          sex = as.matrix(as.numeric(sex == "M"))
          dimnames(sex) = list(rownames(pheno), "sex")
        } # if(length(grep("^sex$", colnames(pheno), ignore.case = TRUE)) > 0)

        x.qtl = qtl.qtlrel(pheno = p[keep,drop = FALSE], probs = mfprobs[keep,,],
                K = K[["X"]][keep,keep], snps = snps[xchr,,drop = FALSE])

      } else {

        keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))

        if(missing(intcovar)) {

          # Additive covariates only.
          x.qtl = fast.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                  K = K[["X"]][keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                  snps = snps[xchr,])

        } else {

          # Additive & interactive covariates.
              x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                    K = K[["X"]][keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                    intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])
        } # else
      } # else

      auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
      auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)

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

