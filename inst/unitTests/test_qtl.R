# QTL mapping unit tests.
test_qtl_no_kinship_no_covar = function() {
  library(MUGAExampleData)
  library(QTLRel)

  # Load in phenotype and genotype data.
  data(pheno)
  data(model.probs)
  probs = model.probs

  # Load in SNPs.
  load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))

  snps = muga_snps[muga_snps[,2] == 1,]
  rm(muga_snps)
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,snps[,1]]

  pheno.name = "WBC1"
  qtl = scanone(pheno = pheno, pheno.col = pheno.name, probs = probs, 
                snps = snps)

  pheno = pheno[!is.na(pheno[,pheno.name]),]
  pheno = pheno[rownames(pheno) %in% rownames(probs),]
  probs = probs[rownames(probs) %in% rownames(pheno),,]
  probs = probs[rownames(pheno),,]
  stopifnot(all(rownames(pheno) == rownames(probs)))

  # Null model.
  null.ss = sum((pheno[,pheno.name] - mean(pheno[,pheno.name]))^2)
  lod = rep(0, nrow(snps))
  pheno.column = which(colnames(pheno) == pheno.name)

  # Full model.
  tmpx = cbind(rep(1, nrow(pheno)), probs[,,1])
  for(s in 1:nrow(snps)) {
    tmpx[,-1] = probs[,,s]
    lod[s] = sum(qr.resid(qr(tmpx), pheno[,pheno.column])^2)
  } # for(s)
  lod = -nrow(pheno) * log10(lod / null.ss) / 2

  checkEqualsNumeric(target = lod, current = qtl$lod$A[,7])

} # test_qtl_no_kinship_no_covar()



test_qtl_no_kinship_with_covar = function() {
  library(MUGAExampleData)
  library(QTLRel)
  
  # Load in phenotype and genotype data.
  data(pheno)
  data(model.probs)
  probs = model.probs

  # Load in SNPs.
  load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))

  snps = muga_snps[muga_snps[,2] == 1,]
  rm(muga_snps)
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,snps[,1]]
  covar = matrix(as.numeric(c(pheno$Sex == "M", pheno$Diet == "hf")),
          ncol = 2, dimnames = list(rownames(pheno), c("sex", "diet")))

  pheno.name = "WBC1"
  qtl = scanone(pheno = pheno, pheno.col = pheno.name, probs = probs, 
                addcovar = covar, snps = snps)

  pheno = pheno[!is.na(pheno[,pheno.name]),]
  pheno = pheno[rownames(pheno) %in% rownames(probs),]
  probs = probs[rownames(probs) %in% rownames(pheno),,]
  probs = probs[rownames(pheno),,]
  covar = covar[rownames(pheno),]
  stopifnot(all(rownames(pheno) == rownames(probs)))
  stopifnot(all(rownames(pheno) == rownames(covar)))

  # Null model.
  null.ss = sum(residuals(lm(pheno[,pheno.name] ~ covar))^2)
  lod = rep(0, nrow(snps))
  pheno.column = which(colnames(pheno) == pheno.name)
  tmpx = cbind(rep(1, nrow(covar)), covar, probs[,,1])

  # Full model.
  for(s in 1:nrow(snps)) {
    tmpx[,-1:-3] = probs[,,s]
    qrx = qr(tmpx)
    lod[s] = sum(qr.resid(qrx, pheno[,pheno.column])^2)
  } # for(s)
  lod = -nrow(pheno) * log10(lod / null.ss) / 2

  checkEqualsNumeric(target = lod, current = qtl$lod$A[,7])

} # test_qtl_no_kinship_with_covar()


test_qtl_with_kinship_no_covar = function() {

  # Load in phenotype and genotype data.
  data(pheno)
  data(model.probs)
  probs = model.probs

  # Load in SNPs.
  load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))

  snps = muga_snps[muga_snps[,2] == 1,]
  rm(muga_snps)
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,snps[,1]]
  K = kinship.probs(probs = probs, snps = snps, bychr = FALSE)

  pheno.name = "WBC1"
  qtl = scanone(pheno = pheno, pheno.col = pheno.name, probs = probs, 
                K = K, snps = snps)

  pheno = pheno[!is.na(pheno[,pheno.name]),]
  pheno = pheno[rownames(pheno) %in% rownames(probs),]
  probs = probs[rownames(probs) %in% rownames(pheno),,]
  probs = probs[rownames(pheno),,]
  K = K[rownames(pheno), rownames(pheno)]
  stopifnot(all(rownames(pheno) == rownames(probs)))
  stopifnot(all(rownames(pheno) == rownames(K)))

  # Estimate variance components and create correction matrix.
  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3],
               snp = snps[,1])
  vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
              EE = diag(nrow(pheno)))
  vc = estVC(y = pheno[,pheno.name], v = vTmp)
  corrMat = vc$par["AA"] * vTmp$AA + vc$par["EE"] * vTmp$EE
  rm(vc, vTmp, prdat)
  eig = eigen(corrMat, symmetric = TRUE)
  if(any(eig$values <= 0)) {
     stop("The covariance matrix is not positive definite")
  } # if(any(eig$values <= 0))
  corrMat = eig$vectors %*% diag(1.0 / sqrt(eig$values)) %*% t(eig$vectors)
  dimnames(corrMat) = list(rownames(pheno), rownames(pheno))

  # Null model.
  ph = corrMat %*% pheno[,pheno.name]
  null.ss = sum(residuals(lm(ph ~ 1))^2)
  lod = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), 9, dimnames = list(snps[,1], c("Intercept", LETTERS[1:8])))
  tmpx = cbind(rep(1, nrow(pheno)), probs[,,1])

  # Full model.
  for(s in 1:nrow(snps)) {
    tmpx[,-1] = probs[,,s]
    tmpxcorr = corrMat %*% tmpx
    qrx = qr(tmpxcorr)
    lod[s] = sum(qr.resid(qrx, ph)^2)
    coef[s,] = qr.coef(qrx, ph)
  } # for(s)
  lod = -nrow(pheno) * log10(lod / null.ss) / 2
  coef[,ncol(coef)] = 0
  coef[,-1] = coef[,-1] + coef[,1]
  coef = coef[,-1]
  coef[,-1] = coef[,-1] - coef[,1]

  checkEqualsNumeric(target = lod, current = qtl$lod$A[,7],
                     tolerance = 1e-7)
  checkEqualsNumeric(target = coef, current = qtl$coef$A)

} # test_qtl_with_kinship_no_covar()



test_qtl_with_kinship_with_covar = function() {
  library(MUGAExampleData)
  library(QTLRel)

  # Load in phenotype and genotype data.
  data(pheno)
  data(model.probs)
  probs = model.probs

  # Load in SNPs.
  load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))

  snps = muga_snps[muga_snps[,2] == 1,]
  rm(muga_snps)
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,snps[,1]]
  covar = matrix(as.numeric(c(pheno$Sex == "M", pheno$Diet == "hf")),
          ncol = 2, dimnames = list(rownames(pheno), c("sex", "diet")))
  K = kinship.probs(probs = probs, snps = snps, bychr = FALSE)

  pheno.name = "WBC1"
  qtl = scanone(pheno = pheno, pheno.col = pheno.name, probs = probs, 
                K = K, addcovar = covar, snps = snps)

  pheno = pheno[!is.na(pheno[,pheno.name]),]
  pheno = pheno[rownames(pheno) %in% rownames(probs),]
  probs = probs[rownames(probs) %in% rownames(pheno),,]
  probs = probs[rownames(pheno),,]
  covar = covar[rownames(pheno),]
  K = K[rownames(pheno), rownames(pheno)]
  stopifnot(all(rownames(pheno) == rownames(probs)))
  stopifnot(all(rownames(pheno) == rownames(K)))
  stopifnot(all(rownames(pheno) == rownames(covar)))

  # Estimate variance components and create correction matrix.
  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3],
               snp = snps[,1])
  vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
              EE = diag(nrow(pheno)))
  vc = estVC(y = pheno[,pheno.name], x = covar, v = vTmp)
  corrMat = vc$par["AA"] * vTmp$AA + vc$par["EE"] * vTmp$EE
  rm(vc, vTmp, prdat)
  eig = eigen(corrMat, symmetric = TRUE)
  if(any(eig$values <= 0)) {
     stop("The covariance matrix is not positive definite")
  } # if(any(eig$values <= 0))
  corrMat = eig$vectors %*% diag(1.0 / sqrt(eig$values)) %*% t(eig$vectors)
  dimnames(corrMat) = list(rownames(pheno), rownames(pheno))

  # Null model.
  ph = corrMat %*% pheno[,pheno.name]
  tmpcovar = corrMat %*% covar
  null.ss = sum(residuals(lsfit(x = tmpcovar, y = ph))^2)
  lod = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), 11, dimnames = list(snps[,1], c("Intercept", 
         "sex", "diet", LETTERS[1:8])))
  tmpx = cbind(rep(1, nrow(pheno)), covar, probs[,,1])
  probs.rng = (2 + ncol(covar)):ncol(tmpx)

  # Full model.
  for(s in 1:nrow(snps)) {
    tmpx[,probs.rng] = probs[,,s]
    qrx = qr(corrMat %*% tmpx)
    lod[s] = sum(qr.resid(qrx, ph)^2)
    coef[s,] = qr.coef(qrx, ph)
  } # for(s)
  lod = -nrow(pheno) * log10(lod / null.ss) / 2
  coef[,ncol(coef)] = 0
  coef[,-1:-3] = coef[,-1:-3] + coef[,1]
  coef[,1] = coef[,colnames(coef) == "A"]
  coef = coef[,colnames(coef) != "A"]
  coef[,-1:-3] = coef[,-1:-3] - coef[,1]

  checkEqualsNumeric(target = lod, current = qtl$lod$A[,7],
                     tolerance = 1e-6)
  checkEqualsNumeric(target = coef, current = qtl$coef$A)

} # test_qtl_with_kinship_with_covar()


