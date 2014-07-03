scanone = function(pheno, pheno.col = 1, probs, K, addcovar, intcovar, snps,
          model = c("additive", "dominance", "full")) {

  model = match.arg(model)

  if(is.null(rownames(pheno))) {
    stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
  } # if(is.null(rownames(pheno)))

  snps[,2] = as.character(snps[,2])

  num.auto = get.num.auto(snps)

  # Convert phenotype names to phenotype column numbers.
  if(is.character(pheno.col)) {
    pheno.col = match(pheno.col, colnames(pheno))
  } # if(is.character(pheno.col))

  # Match up the sample names in the phenotype and genotype matrices.
  pheno = pheno[rownames(pheno) %in% dimnames(probs)[[1]],,drop = FALSE]
  probs = probs[dimnames(probs)[[1]] %in% rownames(pheno),,]
  probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,]

  if(any(dim(probs) == 0)) {
    stop(paste("There are no matching samples in pheno and probs. Please",
         "verify that the sample IDs in rownames(pheno) match the sample",
         "IDs in dimnames(probs)[[1]]."))
  } # if(any(dim(probs) == 0))

#  probs = filter.geno.probs(probs)
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]

  if(!missing(K)) {
     K = K[rownames(K) %in% rownames(pheno), colnames(K) %in% rownames(pheno)]
     K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), colnames(K))]
  } # if(!missing(K))

  if(!missing(addcovar)) {

    if(is.null(rownames(addcovar))) {
      stop("rownames(addcovar) is null. The sample IDs must be in rownames(addcovar).")
    } # if(is.null(rownames(addcovar)))

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

    if(!missing(K)) {
      K = K[rownames(K) %in% rownames(pheno), colnames(K) %in% rownames(pheno)]
      K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), colnames(K))]
    } # if(!missing(K))

    if(is.null(colnames(addcovar))) {
      colnames(addcovar) = paste("addcovar", 1:ncol(addcovar), sep = ".")
    } # if(is.null(colnames(addcovar)))

    addcovar = as.matrix(addcovar)
  } # if(!is.null(addcovar))

  if(!missing(intcovar)) {
    if(is.null(rownames(intcovar))) {
      stop("rownames(intcovar) is null. The sample IDs must be in rownames(intcovar).")
    } # if(is.null(rownames(intcovar)))
    intcovar = as.matrix(intcovar)
    intcovar = intcovar[rownames(intcovar) %in% rownames(pheno),,drop = FALSE]
    intcovar = intcovar[match(rownames(pheno), rownames(intcovar)),,drop = FALSE]
    if(is.null(colnames(intcovar))) {
      colnames(intcovar) = paste("intcovar", 1:ncol(intcovar), sep = ".")
    } # if(is.null(colnames(intcovar)))
  } # if(!is.null(intcovar))

  xchr = which(snps[,2] %in% "X")
  sex = NULL
  if(length(xchr) > 0) {
    sex.col = grep("sex", colnames(pheno), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      if(!missing(addcovar)) {
        sex.col = grep("sex", colnames(addcovar), ignore.case = TRUE)
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
      if(missing(K)) {
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                   snps = snps[auto,])
      } else {
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                   K = K[keep,keep], snps = snps[auto,])
      } # else

    } else {

      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                   rowSums(is.nan(addcovar)) == 0 &
                                   rowSums(is.infinite(addcovar)) == 0))
      if(missing(intcovar)) {

        # Additive covariates only.
        if(missing(K)) {

          auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                     addcovar = addcovar[keep,,drop = FALSE],
                     snps = snps[auto,])
        }  else {

          auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                     K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                     snps = snps[auto,])

        } # else
      } else {

        if(missing(K)) {

          auto.qtl = qtl.qtlrel(pheno = p[keep], probs = probs[keep,,auto],
                     addcovar = addcovar[keep,,drop = FALSE], 
                     intcovar = intcovar[keep,,drop = FALSE], 
                     snps = snps[auto,])

        } else {

          auto.qtl = qtl.qtlrel(pheno = p[keep], probs = probs[keep,,auto],
                     K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                     intcovar = intcovar[keep,,drop = FALSE], snps = snps[auto,])

        } # else
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
        if(length(grep("sex", colnames(pheno), ignore.case = TRUE)) > 0) {
          sex = pheno[,grep("sex", colnames(pheno), ignore.case = TRUE)]
          sex = as.matrix(as.numeric(sex == "M"))
          dimnames(sex) = list(rownames(pheno), "sex")
        } # if(length(grep("sex", colnames(pheno), ignore.case = TRUE)) > 0)

        if(missing(K)) {
          x.qtl = qtl.qtlrel(pheno = p[keep,drop = FALSE], probs = mfprobs[keep,,],
                  addcovar = sex[keep, drop = FALSE], snps = snps[xchr,,drop = FALSE])
        } else {
          x.qtl = qtl.qtlrel(pheno = p[keep,drop = FALSE], probs = mfprobs[keep,,],
                  K = K[keep,keep], addcovar = sex[keep, drop = FALSE], 
                  snps = snps[xchr,,drop = FALSE])
        } # else

      } else {

        keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
        if(missing(intcovar)) {

          # Additive covariates only.
          if(missing(K)) {
            x.qtl = fast.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                    addcovar = addcovar[keep,,drop = FALSE],
                    snps = snps[xchr,])
          } else {
            x.qtl = fast.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                    K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                    snps = snps[xchr,])
          } # else
        } else {
          # Additive & interactive covariates.
          if(missing(K)) {
            x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                  addcovar = addcovar[keep,,drop = FALSE], 
                  intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])
          } else {
            x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                  K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                  intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])
          } # else
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

} # scanone()
