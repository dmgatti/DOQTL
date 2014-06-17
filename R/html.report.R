# Arguments: path: character string containing the directory to write the results to.
#            qtl: list containing the results of scanone (lod & coef).
#            perms: numeric vector containing the permutations for the QTL.
#            merge: boolean, if true, look for corresponding *.Rdata files containing
#                   the names of the qtl in the current working directory and plot
#                   the merge plots. If false (default), do not plot merge analysis.
html.report = function(path, qtl, perms, merge = FALSE) {
  if(missing(qtl)) {
    stop(paste("html.report: qtl argument cannot be missing. Please supply a list",
	     "with each element containing a list with lrs and coef."))
  } # if(missing(qtl))
  sig.qtl = NULL
  # Create the individual QTL pages for each phenotype.
  print("Creating Individual QTL pages...")
  for(i in 1:length(qtl)) {
    print(names(qtl)[i])
    tmp = create.html.page(path = path, qtl = qtl[[i]], pheno.name = names(qtl)[i],
          perms = perms)
    if(!is.null(tmp)) {
      sig.qtl = rbind(sig.qtl, tmp)
    } else {
      sig.qtl = rbind(sig.qtl, c(names(qtl)[i], "No", "signficant", "QTL", rep("", 6)))
    } # else
  } # for(i)
  # Make the main page with a summary table.
  print("Creating main QTL page...")
  p = openPage("index.html")
  hwrite("QTL Summary", p, br = TRUE, center = TRUE, heading = 1)
  sig.qtl = sig.qtl[order(sig.qtl[,1]),]
  write.csv(sig.qtl, "sig.qtl.csv")
  hwrite("QTL Summary File", p, br = TRUE, center = TRUE, link = "sig.qtl.csv")
  tmp = make.names(sig.qtl[,1])
  for(i in 1:nrow(sig.qtl)) {
    sig.qtl[i,1] = hwrite(sig.qtl[i,1], link = paste(tmp[i], "/", tmp[i],
                   ".QTL.html", sep = ""))
  } # for(i)
  hwrite(sig.qtl, p, br = TRUE, row.bgcolor = "#aaaaaa")
  closePage(p)
} # html.report()
##########################################################
# Function to create a QTL results page for one phenotype.
# Returns a data.frame with the significant QTL.
create.html.page = function(path, qtl, pheno.name, perms) {
  curr.name = make.names(pheno.name)
  dir.create(paste(path, "/", curr.name, sep = ""))
  curr.path = paste(path, "/", curr.name, "/", sep = "")
  thr = NULL
  if(!missing(perms)) {
    thr = quantile(perms, c(0.37, 0.9, 0.95))
  } # if(!missing(perms))
    
  p = openPage(paste(curr.path, curr.name, ".QTL.html", sep = ""))
  hwrite(paste(pheno.name, "QTL"), p, br = TRUE, heading = 1)
  sig.qtl = NULL
  if(!is.null(thr)) {
    image.file = paste(curr.name, "QTL.png", sep = "")
    png(paste(curr.path, image.file, sep = ""), width = 1000, height = 800,
        res = 128)
    plot.doqtl(qtl, main = pheno.name, sig.thr = thr,
       sig.col = c("goldenrod", "orange", "red"))
    dev.off()
    hwriteImage(image.file, p, br = TRUE, link = image.file)
    # Significant QTL for this phenotype. We can only determine significant
    # peaks if the user has provided thresholds.
    cn = c("Phenotype", "SNP.ID", "Chr", "Pos (Mb)", "Pos (cM)", "LOD",
           "p-value", "% Var",  "Proximal (Mb)", "Distal (Mb)")
    sig.tmp = rep("", length(cn))
    sig.A.chr = unique(qtl$lod$A[which(qtl$lod$A$lod > thr[1]),2])
    sig.X.chr = unique(qtl$lod$X[which(qtl$lod$X$lod > thr[1]),2])
    if(length(sig.A.chr) > 0 | length(sig.X.chr) > 0) {
      coef.image.files = NULL
      merge.image.files = NULL
      if(length(sig.A.chr) > 0) {
        coef.image.files = paste(curr.name, ".chr", sig.A.chr, ".coef.png",
                           sep = "")
        for(c in 1:length(sig.A.chr)) {
          png(paste(curr.path, coef.image.files[c], sep = ""), width = 1000, 
              height = 800, res = 128)
          coefplot(qtl, chr = sig.A.chr[c], main = pheno.name)
          dev.off()
          interval = bayesint(qtl$lod$A, chr = sig.A.chr[c], expandtomarkers = TRUE)
          sig.tmp[1] = pheno.name
          sig.tmp[2:6] = interval[2,c(1:4,7)]
          sig.tmp[7] = mean(perms >= as.numeric(interval[2,7]))
          sig.tmp[8] = interval[2,5]
          sig.tmp[9:10] = interval[c(1,3),3]
          sig.qtl = rbind(sig.qtl, unlist(sig.tmp))
        } # for(c)
      } # if(length(sig.A.chr) > 0)
      if(length(sig.X.chr) > 0) {
        
        # If we have 7 female coeffeficients, then plot them.
        # TBD: This is DO specific!
        if(all(paste("F", LETTERS[2:8], sep = ".") %in% colnames(qtl$coef$X))) {
          coef.image.files = c(coef.image.files, paste(curr.name, 
                             ".chrX.female.coef.png", sep = ""))
          png(paste(curr.path, coef.image.files[length(coef.image.files)], 
              sep = ""), width = 1000, height = 800, res = 128)
          coefplot(qtl, chr = "X", main = paste(pheno.name, "Females"), sex = "F")
          dev.off()
          interval = bayesint(qtl$lod$X, chr = "X", expandtomarkers = TRUE)
          sig.tmp[1] = pheno.name
          sig.tmp[2:6] = interval[2,c(1:4,7)]
          sig.tmp[7] = mean(perms >= as.numeric(interval[2,7]))
          sig.tmp[8] = interval[2,5]
            sig.tmp[9:10] = interval[c(1,3),3]
          sig.qtl = rbind(sig.qtl, unlist(sig.tmp))
        } # if(all(paste("F", LETTERS[2:8], sep = ".") %in% colnames(qtl$coef$X)))
        # If we have 7 male coefficients, then plot them.
        # TBD: This is DO specific!
        if(all(paste("M", LETTERS[2:8], sep = ".") %in% colnames(qtl$coef$X))) {
          coef.image.files = c(coef.image.files, paste(curr.name, 
                             ".chrX.male.coef.png", sep = ""))
          png(paste(curr.path, coef.image.files[length(coef.image.files)], 
              sep = ""), width = 1000, height = 800, res = 128)
          coefplot(qtl, chr = "X", main = paste(pheno.name, "Males"), sex = "M")
          dev.off()
          interval = bayesint(qtl$lod$X, chr = "X", expandtomarkers = TRUE)
          sig.tmp[1] = pheno.name
          sig.tmp[2:6] = interval[2,c(1:4,7)]
          sig.tmp[7] = mean(perms >= as.numeric(interval[2,7]))
          sig.tmp[8] = interval[2,5]
            sig.tmp[9:10] = interval[c(1,3),3]
          sig.qtl = rbind(sig.qtl, unlist(sig.tmp))
        } # if(all(paste("M", LETTERS[2:8], sep = ".") %in% colnames(qtl$coef$X)))
        colnames(sig.qtl) = cn
        sig.qtl[,6] = format(as.numeric(sig.qtl[,6]), digits = 4)
        sig.qtl[,8] = format(as.numeric(sig.qtl[,8]) * 100, digits = 4)
        hwrite(sig.qtl, p, br = TRUE, row.bgcolor = "#aaaaaa")
      } # if(length(sig.X.chr) > 0)
      for(i in 1:length(coef.image.files)) {
        hwriteImage(coef.image.files[i], p, br = TRUE, link = coef.image.files[i])
      } # for(c)
    } #     if(length(sig.A.chr) > 0 | length(sig.X.chr) > 0)
  } else {
    plot.doqtl(qtl$lod, main = pheno.name)
  } # else
  closePage(p)
  return(sig.qtl)
} # create.html.page()
