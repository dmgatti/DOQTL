scanone <- 
function(pheno, pheno.col = 1, probs, K, addcovar, intcovar, snps,
         model = c("additive", "dominance", "full")) {

  model = match.arg(model)

  if(is.null(rownames(pheno))) {
    stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
  } # if(is.null(rownames(pheno)))

  if(is.list(K)) {
    chr = unique(snps[,2])
    if(all(!names(K) %in% chr)) {
       stop(paste("All of the chromosomes in the K list are not in snps .",
            "Please supply a list of kinship matrices names for each of",
            "the chromosomes in snps."))
    } # if(all(!chr %in% names(K))
  } # if(is.list(K))

  num.auto = get.num.auto(snps)

  if(is.character(pheno.col)) {
    pheno.col = match(pheno.col, colnames(pheno))
  } # if(is.character(pheno.col))

  # Match up the sample names in the phenotype and genotype matrices.
  pheno = pheno[rownames(pheno) %in% dimnames(probs)[[1]],,drop = F]
  probs = probs[dimnames(probs)[[1]] %in% rownames(pheno),,]
  probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,]

  if(any(dim(probs) == 0)) {
    stop(paste("There are no matching samples in pheno and probs. Please",
         "verify that the sample IDs in rownames(pheno) match the sample",
         "IDs in dimnames(probs)[[1]]."))
  } # if(any(dim(probs) == 0))

  probs = filter.geno.probs(probs)
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
  if(is.list(K)) {
    K = lapply(K, function(z) {
          z = z[rownames(z) %in% rownames(pheno), colnames(z) %in% rownames(pheno)]
          z[match(rownames(pheno), rownames(z)), match(rownames(pheno), colnames(z))]
        }) # lapply()
  } else {
    K = K[rownames(K) %in% rownames(pheno), colnames(K) %in% rownames(pheno)]
    K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), colnames(K))]
  } # else

  xchr = which(snps[,2] %in% "X")
  sex = NULL
  if(length(xchr) > 0) {
    sex.col = grep("sex", colnames(pheno), ignore.case = T)
    if(length(sex.col) == 0) {
      if(!missing(addcovar)) {
        sex.col = grep("sex", colnames(addcovar), ignore.case = T)
      }
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
  } # if(length(xchr) > 0)

  if(!missing(addcovar)) {
    if(is.null(rownames(addcovar))) {
      stop("rownames(addcovar) is null. The sample IDs must be in rownames(addcovar).")
    } # if(is.null(rownames(addcovar)))
    addcovar = data.frame(addcovar, stringsAsFactors = T)
    rn = rownames(addcovar)
    addcovar = lapply(addcovar, as.numeric)
    cn = names(addcovar)
    addcovar = matrix(unlist(addcovar), length(addcovar[[1]]),
               length(addcovar), dimnames = list(rn, cn))
    rownames(addcovar) = rn
    addcovar = addcovar[rownames(addcovar) %in% rownames(pheno),,drop = F]
    addcovar = addcovar[match(rownames(pheno), rownames(addcovar)),,drop = F]

    if(sum(rownames(pheno) %in% rownames(addcovar)) == 0) {
      stop(paste("rownames(pheno) does not contain any sample IDs in",
           "common with rownames(addcovar). Please make sure that the",
           "rownames in pheno and addcovar match."))
    } # if(sum(rownames(pheno) %in% rownames(addcovar)) == 0)

    pheno = pheno[rownames(pheno) %in% rownames(addcovar),,drop = F]
    addcovar = addcovar[rownames(addcovar) %in% rownames(pheno),,drop = F]
    addcovar = addcovar[match(rownames(pheno), rownames(addcovar)),,drop = F]
    probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]
    if(is.list(K)) {
      K = lapply(K, function(z) {
            z[rownames(z) %in% rownames(pheno), colnames(z) %in% rownames(pheno)]
          }) # lapply()
    } else {
      K = K[rownames(K) %in% rownames(pheno), colnames(K) %in% rownames(pheno)]
    } # else

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
    intcovar = intcovar[rownames(intcovar) %in% rownames(pheno),,drop = F]
    intcovar = intcovar[match(rownames(pheno), rownames(intcovar)),,drop = F]
    if(is.null(colnames(intcovar))) {
      colnames(intcovar) = paste("intcovar", 1:ncol(intcovar), sep = ".")
    } # if(is.null(colnames(intcovar)))
  } # if(!is.null(intcovar))

  # Make the results list.
  retval = as.list(1:length(pheno.col))
  names(retval) = colnames(pheno)[pheno.col]
  index = 1

  for(i in pheno.col) {
    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)

    keep = which(!is.na(p))

    # Autosomes
    auto = which(snps[,2] %in% 1:num.auto)

    if(missing(addcovar)) {
      # No covariates.
      if(is.list(K)) {
        for(i in 1:length(K)) {
          auto.qtl = fast.qtlrel(pheno = p[keep], prob = probs[keep,,auto], 
                     K = K[keep,keep], snps = snps[auto,])        
        } # for(i)
      } else {
        auto.qtl = fast.qtlrel(pheno = p[keep], prob = probs[keep,,auto], 
                   K = K[keep,keep], snps = snps[auto,])
      } # else
    } else {
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
      if(missing(intcovar)) {
        # Additive covariates only.
        auto.qtl = fast.qtlrel(pheno = p[keep], prob = probs[keep,,auto], 
                   K = K[keep,keep], addcovar = addcovar[keep,,drop = F],
                   snps = snps[auto,])
      } else {
        auto.qtl = qtl.qtlrel(pheno = p[keep], prob = probs[keep,,auto],
                   K = K[keep,keep], addcovar = addcovar[keep,,drop = F], 
                   intcovar = intcovar[keep,,drop = F], snps = snps[auto,], 
                   test = "None")
      } # else
    } # else
    auto.qtl = list(lod = list(A = auto.qtl$lod), coef = list(A = auto.qtl$coef))

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
        mfprobs = probs
        dimnames(mfprobs)[[2]] = paste("F", dimnames(mfprobs)[[2]], sep = ".")
      } else if(length(males) > 0) {
        mfprobs = probs
        dimnames(mfprobs)[[2]] = paste("M", dimnames(mfprobs)[[2]], sep = ".")
      } # else if(length(males) > 0)

      if(missing(addcovar)) {
        # No covariates, but we need sex for the null hypothesis.
        if(length(grep("sex", colnames(pheno), ignore.case = T)) > 0) {
          sex = pheno[,grep("sex", colnames(pheno), ignore.case = T)]
          sex = as.matrix(sex)
          dimnames(sex) = list(rownames(pheno), "sex")
        } # if(length(grep("sex", colnames(pheno), ignore.case = T)) > 0)
        x.qtl = qtl.qtlrel(pheno = p[keep,drop = F], prob = mfprobs[keep,,],
                K = K[keep,keep], addcovar = sex[keep, drop = F], 
                snps = snps[xchr,,drop = F], test = "None")
      } else {
        keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
        if(missing(intcovar)) {
          # Additive covariates only.
          x.qtl = fast.qtlrel(pheno = p[keep], prob = mfprobs[keep,,], 
                  K = K[keep,keep], addcovar = addcovar[keep,,drop = F],
                  snps = snps[xchr,])
        } else {
          x.qtl = qtl.qtlrel(pheno = p[keep], prob = mfprobs[keep,,],
                K = K[keep,keep], addcovar = addcovar[keep,,drop = F], 
                intcovar = intcovar[keep,,drop = F], snps = snps[xchr,], 
                test = "None")
        } # else
      } # else

      auto.qtl$lod = list(A = auto.qtl$lod$A, X = x.qtl$lod)
      auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)
    } # if(length(xchr) > 0)

    retval[[index]] = auto.qtl
    index = index + 1
  } # for(i)

  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)

  class(retval) = c(class(retval), "DOQTL")
  attr(retval, "model") = "additive"

  return(retval)

} # scanone()

