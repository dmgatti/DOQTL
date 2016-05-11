################################################################################
# Main function for running the calculating genotype probabilities.
# The package contains routines to estimate recombination frequencies for the
# Collaborative Cross and Diversity Outbred mice. Other mouse populations or
# organisms will need to supply snps, founder data and transition probabilities.
# Arguments: data: list containing the data needed to map in these samples.
#            chr: character vector, with chromosomes to run.  Must match
#                 the chromosome IDs in the snps table.  all means
#                 run all.
#            output.dir: character, with the file directory for the final
#                        theta & rho mean and variance files.
#            plot: boolean that is true if the user would like to plot the
#                  genotypes after reconstruction.
#            array: character, with the type of array used to genotype the 
#                   samples. The muga and megamuga options are for the Mouse
#                   Universal Genotyping Arrays. Options are 'muga', 'megamuga',
#                   or 'other'.
#            sampletype: character, with the type of samples. Options are 'DO',
#                        'CC', 'DOF1', 'HS' or 'other'.
#            method: character, that inidicates the method of genotype reconstruction.
#                    'intensity' uses X and Y array intensities and 'allele' uses
#                    the genotype calls.
#            founders: list, containing founder data. Must contain allele calls
#                         or X & Y intensites, founder names, codes and sexes.
#            transprobs: list, containing transition probabilities between each pair
#                        of SNPs on each chromosome.
#            snps: data.frame with SNP IDs, chr, Mb and cM positions. For custom
#                  arrays.
calc.genoprob = function(data, chr = "all", output.dir = ".", plot = TRUE, 
                         array = c("gigamuga", "megamuga", "muga", "other"),
                         sampletype = c("DO", "CC", "DOF1", "HS", "HSrat", "other"),
                         method = c("intensity", "allele"), clust = c("mclust", "pamk"),
                         founders, transprobs, snps) {

  if(missing(data)) {
    stop(paste("data is missing. Please supply data to process."))
  } # if(missing(data))

  if(is.null(data$sex)) {
    stop("data$sex cannot be null, even if all animals are the same sex.")
  } else {
    if(!is.vector(data$sex)) { 
      stop("data$sex must be a vector containing only M and F with sample IDs in names.")
    } # if(!is.vector(data$sex))

    # Remove samples with no sex.
    if(any(is.na(data$sex))) {
      wh = which(is.na(data$sex))

      warning(paste("Some samples have sex = NA. These will be removed.", 
              names(sex)[wh]))

      sex = sex[-wh]
      gen = gen[-wh]
      if("x" %in% names(data)) {
        x = x[-wh,]
        y = y[-wh,]
      } else {
        geno = geno[-wh,]
      } # else
    } # if(any(is.na(data$sex))
  } # else

  if(!is.null(data$gen)) {
    if(!is.vector(data$gen)) { 
      stop("data$gen must be a vector containing only numbers with sample IDs in names.")
    } # if(!is.vector(data$gen))
  } # if(!is.null(data$gen))

  method = match.arg(method)
  sampletype = match.arg(sampletype)
  array = match.arg(array)
  clust = match.arg(clust)
  attr(data, "method") = method
  attr(data, "array")  = array
  attr(data, "sampletype") = sampletype

  if(method == "allele") {

    if(any(!c("geno", "sex", "gen") %in% names(data))) {
      stop(paste("\'data\' does not contain an element called",
           c("geno", "sex", "gen")[!c("geno", "sex", "gen") %in% names(data)],
           ".\nPlease make sure that \'data\' contains three elements:",
           "\'geno\', \'sex\' and \'gen\'."))
    } # if(any(!c("geno", ...

  } else if(method == "intensity") {
    if(any(!c("x", "y", "sex", "gen") %in% names(data))) {
      stop(paste("\'data\' does not contain an element called",
           c("x", "y", "sex", "gen")[!c("x", "y", "sex", "gen") %in% names(data)],
           ".\nPlease make sure that \'data\' contains three elements:",
           "\'x\', \'y\', \'sex\' and \'gen\'."))
    } # if(any(!c("x", "y", ...
  } # else

  # If the sampletype is CC or DO and array is muga or megamuga, then get the 
  # founder data.
############################
############################
###       CC or DO       ###
############################
############################
  if(sampletype == "CC" | sampletype == "DO") {

    ### MUGA Series Arrays ###
    if(array %in% c("gigamuga", "megamuga", "muga")) { 

      founders = read.muga.data(array = array, method = method)
      snps = founders$snps

      if(array == "muga") {

        # Get the MUGA SNPs and subset them to keep the ones that we use for 
        # genotyping.
        snps = snps[snps[,1] %in% muga.snps.to.keep,]
        snps = snps[!is.na(snps[,4]),1:4]

      ### MegaMUGA ###
      } else if(array == "megamuga") {

        # Get the MUGA SNPs and subset them to keep the ones that we use for 
        # genotyping.
        snps = snps[snps$tier <= 2,]
        snps = snps[!is.na(snps[,4]),1:4]

        # Remove founders with too much missing data.
        if(method == "allele") {

          # Remove founders with too much missing data.
          remove = which(rowMeans(founders$geno == "N") > 0.05)
          founders$geno = founders$geno[-remove,]
          founders$sex  = founders$sex[-remove]
          founders$code = founders$code[-remove]

        } else if(method == "intensity") {

          # Remove founders with too much missing data.
          remove = which(rowSums(founders$x < 0) > 1000)
          if(length(remove) > 0) {
            founders$x = founders$x[-remove,]
            founders$y = founders$y[-remove,]
            founders$sex  = founders$sex[-remove]
            founders$code = founders$code[-remove]
          } # if(length(remove) > 0)

        } # else

      ### GigaMUGA ###
      } else if(array == "gigamuga") {

        snps = snps[snps$tier <= 2 & snps$is.biallelic == TRUE,]

        # Make sure that the marker positions always increase.
        for(i in 2:nrow(snps)) {
          if(snps[i-1,2] == snps[i,2]) {  
            if(snps[i,3] <= snps[i-1,3]) {
              snps[i,3] = snps[i-1,3] + 0.001
            } # if(snps[i,3] <= snps[i-1,3])
          } # if(snps[i-1,2] == snps[i,2])
        } # for(i)

      } # else if(array == "gigamuga")

      founders = founders[names(founders) != "snps"]

      # Keep only the DO founders.
      founders = keep.do.founders(founders)

      # Add the genotype states.
      states = create.genotype.states(LETTERS[1:8], sampletype)

      if(method == "allele") {
        founders = list(geno = founders$geno, sex = founders$sex, code = founders$code,
                   states = states)
      } else {
        founders = list(x = founders$x, y = founders$y, sex = founders$sex,
                   code = founders$code, states = states)
      } # else
      attr(founders, "method") = method

    ### Custom Array ###
    } else {

      # array = 'other'
      # We have a custom array. The user must supply SNPs, founder data.
      if(missing(founders)) {
        stop(paste("calc.genoprob: founders is missing. You must supply",
             "founder genotypes when using a custom array."))
      } # if(missing(founders))

      if(missing(snps)) {
        stop(paste("calc.genoprob: snps is missing. You must supply",
             "array SNPs when using a custom array."))
      } # if(missing(founders))

      ### Allele Call ###
      if(method == "allele") {

        founders = list(geno = founders$geno, sex = founders$sex,
                   code = founders$code, states = 
                   create.genotype.states(LETTERS[1:8]))

        # This is needed because the as.matrix() function keeps stripping the
        # rownames!
        rn = rownames(data$geno)
        data$geno = as.matrix(data$geno)
        rownames(data$geno) = rn

        if(!all(snps[,1] == colnames(data$geno))) {
          stop(paste("The SNP names in the snps argument do not match the SNP",
               "names in the data."))
        } # if(!all(snps[,1] == colnames(data$geno)))

        if(!all(snps[,1] == colnames(founders$geno))) {
          stop(paste("The SNP names in the snps argument do not match the SNP",
               "names in the founders."))
        } # if(!all(snps[,1] == colnames(founders$geno)))

      ### Intensity ###

      } else if(method == "intensity") {

        founders = list(x = founders$x, y = founders$y,
                   code = founders$code, sex = founders$sex, states = 
                   create.genotype.states(LETTERS[1:8]))

        if(!all(snps[,1] == colnames(data$x))) {
          stop(paste("The SNP names in the snps argument do not match the SNP",
               "names in the data."))
        } # if(!all(snps[,1] == colnames(data$x)))

        if(!all(snps[,1] == colnames(founders$x))) {
          stop(paste("The SNP names in the snps argument do not match the SNP",
               "names in the founders."))
        } # if(!all(snps[,1] == colnames(founders$x)))

      } # else if(method == "intensity")
    } # else

    # Get the outbreeding generation for the DO.
    if(sampletype == "DO") {
      if(!"gen" %in% names(data)) {
        stop(paste("The data list does not contain an element called 'gen'.",
             "Please add the DO outbreeding generation of each sample when",
             "sampletype = 'DO'."))
      } else {
        # Convert the DO# values into just numbers.
        nm = names(data$gen)
        data$gen = as.numeric(sub("^DO", "", data$gen))
        names(data$gen)= nm
      } # else
    } # if(sampletype == "DO")

###############
###  DO/F1  ###
###############
  } else if (sampletype == "DOF1") {

    # Users must provide genotypes for the mutant founders.
    mut.founders = founders

    # Get the DO founders from the website.
    founders = read.muga.data(array = array, method = method)
    snps = founders$snps

    ### MegaMUGA ###
    if(array == "megamuga") {

      # Keep only the collaborative cross SNPs.
#      snps = snps[snps$Collaborative.Cross == 1 | snps$MUGA == 1 |
#	              snps$C57BL.6 == 1,]
      snps = snps[snps$Collaborative.Cross == 1 | snps$Chr.Y == 1,]
      snps = snps[!is.na(snps[,4]) | snps[,2] == "Y",]
      snps = snps[,1:4]

    } else if (array == "gigamuga") {

      snps = snps[grep("^(B6|JAX|ICR|UNC|Xi)", snps[,1]),]

      # Make sure that the marker positions always increase.
      for(i in 2:nrow(snps)) {
        if(snps[i-1,2] == snps[i,2]) {  
          if(snps[i,3] <= snps[i-1,3]) {
            snps[i,3] = snps[i-1,3] + 0.001
          } # if(snps[i,3] <= snps[i-1,3])
        } # if(snps[i-1,2] == snps[i,2])
      } # for(i)

    } else {

      stop(paste("calc.genoprob: Unsupported array platform:", array))

    } # else

      ### Allele Call ###
      if(method == "allele") {

        # Remove founders with too much missing data.
        remove = which(rowSums(founders$geno == "N") > 5000)
        founders$geno = founders$geno[-remove,]
        founders$sex  = founders$sex[-remove]
        founders$code = founders$code[-remove]

        founders = list(geno = t(founders$geno[rownames(founders$geno) %in% snps[,1],]),
                        sex = founders$sex, code = founders$code,
                        states = create.genotype.states(LETTERS[1:8]),
                        direction = founders$direction)
        keep = which(!is.na(founders$code))
        founders$geno = founders$geno[keep,]
        founders$sex  = founders$sex[keep]
        founders$code = founders$code[keep]

        # Add the mutant founders to the founders.
        founders$geno = cbind(mut.founders$geno, founders$geno)

      ### Intensity ###
      } else if(method == "intensity") {

        # Remove founders with too much missing data.
        remove = which(rowSums(founders$x < 0) > 1000)
        founders$x = founders$x[-remove,]
        founders$y = founders$y[-remove,]
        founders$sex  = founders$sex[-remove]
        founders$code = founders$code[-remove]

        # Remove founders with NA founder codes (i.e. non-DO founders)
        remove = which(is.na(founders$code))
        founders$x = founders$x[-remove,]
        founders$y = founders$y[-remove,]
        founders$sex  = founders$sex[-remove]
        founders$code = founders$code[-remove]

        # Keep only the inbred founders and discard the F1s.
        founders$code = founders$code[founders$code %in% paste(LETTERS[1:8], LETTERS[1:8], sep = "")]
        founders$x = founders$x[rownames(founders$x) %in% names(founders$code),]
        founders$y = founders$y[rownames(founders$y) %in% names(founders$code),]
        founders$sex = founders$sex[names(founders$sex) %in% names(founders$code)]

        stopifnot(all(colnames(founders$x) == colnames(founders$x)))
        stopifnot(all(colnames(founders$y) == colnames(founders$y)))

        # Add the mutant founders to the founders.
        founders$geno = cbind(mut.founders$geno, founders$geno)

        if(!all(colnames(data$x) == colnames(founders$x))) {
          stop("calc.genoprob: SNP names in data$x and founders$x do not match for MegaMUGA.")
        } # if(!all(colnames(data$x) == colnames(founders$x)))

        if(!all(colnames(data$y) == colnames(founders$y))) {
          stop("calc.genoprob: SNP names in data$y and founders$y do not match for MegaMUGA.")
        } # if(!all(colnames(data$y) == colnames(founders$y)))

      } # else

      founders$sex  = cbind(mut.founders$sex, founders$sex)
      founders$code = cbind(mut.founders$code, founders$code)

      states = paste(LETTERS[1:8], "I", sep = "")
      if(founders$direction == "DOxMUT") {
        states = list(auto = states, X = list(F = states, M = LETTERS[1:8]),
                 founders = LETTERS[1:9])
      } else if(founders$direction == "MUTxDO") {
        states = list(auto = states, X = list(F = states, M = "I"),
                 founders = LETTERS[1:9])
      } else {
        stop(paste("calc.genoprob: founders$direction =", founders$direction, 
             ". founders$direction must be either MUTxDO or DOxMUT. See",
             "help(calc.genoprob) for more information."))
      } # else

      founders = list(x = x, y = y, sex = sex, code = code,
                      states = states, direction = founders$direction)
      attr(founders, "method") = method

  } else if(sampletype == "HS") {

    # In this case, we require the user to supply founders and markers.

  } else {
###############################
### All other sample types. ###
###############################
    # We require the user to provide the founder data and markers.
    if(missing(founders)) {
      stop("Founder data missing. Founder data must be provided for non-DO samples.")
    } # if(missing(founders))

    if(missing(snps)) {
      stop("Markes missing. Markers must be provided for non-DO samples.")
    } # if(missing(snps))

    if(method == "allele") {
      fnames = c("geno", "sex", "code", "states")
      if(any(!names(founders) %in% fnames)) {
        wh = fnames[!fnames %in% names(founders)]
        print(paste("Required list elements called", wh,
              "were not found in the founder data. Please include these elements."))
      } # if(any(!names(founders) %in% fnames))
    } else {
      fnames = c("x", "y", "sex", "code", "states")
      if(any(!names(founders) %in% fnames)) {
        wh = fnames[!fnames %in% names(founders)]
        print(paste("Required list elements called", wh,
              "were not found in the founder data. Please include these elements."))
      } # if(any(!names(founders) %in% fnames))
    } # else

  } # else

  # Set sex identifiers to upper case.
  data$sex = toupper(data$sex)
  founders$sex = toupper(founders$sex)

  # Fill in the SNP cM values for which consecutive SNPs have the same cM value,
  # even though they are at different Mb positions.
#  snps = fill.in.snps(snps)

  # Fill in any missing F1s.
  attr(founders, "method") = attr(data, "method")
  founders = add.missing.F1s(founders, snps, sampletype)

  # Add a slash to the output directory, if required.
  output.dir = add.slash(output.dir)

  # Set the transition probability function.
  trans.prob.fxn = do.trans.probs

  if(sampletype == "CC") {
    trans.prob.fxn = cc.trans.probs
  } else if(sampletype == "DOF1") {
    trans.prob.fxn = dof1.trans.probs
  } else if(sampletype == "HS" | sampletype == "HSrat") {
    trans.prob.fxn = hs.trans.probs
  } else if(sampletype == "other") {
    trans.prob.fxn = generic.trans.probs
  } # else

  chr = as.character(chr)

  # Verify that SNPs are sorted.
  tmp = split(snps[,1:4], snps[,2])
  tmp = lapply(tmp, function(z) { diff(z[,4]) })
  if(any(unlist(sapply(tmp, "<", 0)), na.rm = TRUE)) {
    wh = sapply(tmp, "<", 0)
    print(paste("The SNPs are not sorted in increasing order on Chr",
          paste(names(wh)[which(sapply(wh, sum) > 0)]), ". Please",
          "sort the SNPs on each chromosome such that the cM value",
          "increases monotonically."))
  } # if(any(sapply(tmp, "<" 0))

  # If chr == all, then change the chr value to the chromosomes we have.
  if(chr == "all") {
    chr = unique(snps[,2])
  } # if(chr == "all")

  # Subset the markers to include only the chromosomes that we are running.
  snps = snps[snps[,2] %in% chr,]

  # Synchronize the markers in the allele calls and snps.
  tmp = synchronize.snps(snps, data, founders)
  data = tmp$data
  founders = tmp$founders
  snps = tmp$markers
  rm(tmp)
  gc()

  # Make the names unique.
  names(data$sex) = make.unique(make.names(names(data$sex)))
  if(attr(data, "sampletype") %in% c("DO", "DOF1", "HS", "HSrat")) {
    names(data$gen) = make.names(names(data$gen))
  } # if(attr(data, "sampletype") == "DO")

  if(method == "allele") {

     # Verify that the allele calls are the same in the data and
     # founders.
     data$geno[data$geno == "-"] = "N"
     check.alleles(data = data, founders = founders)

     # Convert the sample names into names that can be written out 
     # as files.
     rownames(data$geno) = make.unique(make.names(rownames(data$geno)))

     if(attr(data, "sampletype") == "DO") {
       names(data$gen) = make.names(names(data$gen))
     } # if(attr(data, "sampletype") == "DO")

     # Split up the data by chromosome.
     tmpdir = tempdir()
     data$geno = split(data.frame(t(data$geno)), snps[,2])
     founders$geno = split(data.frame(t(founders$geno)), snps[,2])
     data$geno = lapply(data$geno, function(z) { as.matrix(t(z)) })
     founders$geno = lapply(founders$geno, function(z) { as.matrix(t(z)) })

     # Write data out to temporary files.
     savefxn = function(d, f) { save(d, file = f) }
     mapply(savefxn, data$geno, paste(tmpdir, "/data_geno_chr", names(data$geno),
            ".Rdata", sep = ""))
     mapply(savefxn, founders$geno, paste(tmpdir, "/founders_geno_chr", 
            names(founders$geno), ".Rdata", sep = ""))
     # Replace the data with temporary file names.
     data$geno = dir(path = tmpdir, pattern = "data_geno_chr", full.names = TRUE)
     founders$geno = dir(path = tmpdir, pattern = "founders_geno_chr",
                     full.names = TRUE)

     # Split up the SNPs and order them numerically.
     snps.spl = split(snps, snps[,2])
     old.warn = options("warn")$warn
     options(warn = -1)
     snps.spl = snps.spl[order(as.numeric(names(snps.spl)))]
     options(warn = old.warn)
     gc()

     calc.genoprob.alleles(data = data, chr = chr, founders = founders, snps = snps.spl,
                           output.dir = output.dir, trans.prob.fxn = trans.prob.fxn,
                           plot = plot)

  } else if(method == "intensity") {

     if(any(!is.nan(founders$x))) {
       print("NaN values in founder X intensities.")
     } # if(any(!is.nan(founders$x)))

     if(any(!is.nan(founders$y))) {
       print("NaN values in founder Y intensities.")
     } # if(any(!is.nan(founders$y)))

     stopifnot(ncol(data$x) > 0 & ncol(data$y) > 0)
     stopifnot(ncol(founders$x) > 0 & ncol(founders$y) > 0)

     # Convert the sample names into names that can be written out 
     # as files.
     rownames(data$x) = make.unique(make.names(rownames(data$x)))
     rownames(data$y) = make.unique(make.names(rownames(data$y)))

     # Split up the data by chromosome.
     tmpdir = tempdir()
     data$x = split(data.frame(t(data$x)), snps[,2])
     data$y = split(data.frame(t(data$y)), snps[,2])
     founders$x = split(data.frame(t(founders$x)), snps[,2])
     founders$y = split(data.frame(t(founders$y)), snps[,2])
     data$x = lapply(data$x, function(z) { as.matrix(t(z)) })
     data$y = lapply(data$y, function(z) { as.matrix(t(z)) })
     founders$x = lapply(founders$x, function(z) { as.matrix(t(z)) })
     founders$y = lapply(founders$y, function(z) { as.matrix(t(z)) })

     # Write data out to temporary files.
     savefxn = function(d, f) { save(d, file = f) }
     mapply(savefxn, data$x, paste(tmpdir, "/data_x_chr", names(data$x),
            ".Rdata", sep = ""))
     mapply(savefxn, data$y, paste(tmpdir, "/data_y_chr", names(data$y),
            ".Rdata", sep = ""))
     mapply(savefxn, founders$x, paste(tmpdir, "/founders_x_chr", 
            names(founders$x), ".Rdata", sep = ""))
     mapply(savefxn, founders$y, paste(tmpdir, "/founders_y_chr", 
            names(founders$y), ".Rdata", sep = ""))

     # Replace the data with temporary file names.
     data$x = dir(path = tmpdir, pattern = "data_x_chr", full.names = TRUE)
     data$y = dir(path = tmpdir, pattern = "data_y_chr", full.names = TRUE)
     founders$x = dir(path = tmpdir, pattern = "founders_x_chr", full.names = TRUE)
     founders$y = dir(path = tmpdir, pattern = "founders_y_chr", full.names = TRUE)    

     # Split up the SNPs and order them numerically.
     snps.spl = split(snps, snps[,2])
     old.warn = options("warn")$warn
     options(warn = -1)
     snps.spl = snps.spl[order(as.numeric(names(snps.spl)))]
     options(warn = old.warn)
     gc()

     calc.genoprob.intensity(data = data, chr = chr, founders = founders, 
                             snps = snps.spl, output.dir = output.dir, 
                             trans.prob.fxn = trans.prob.fxn, plot = plot, clust = clust)

  } # else if(method = "intensity")

  # Convert the *.genotype.probs.txt files to *.Rdata files.
  create.Rdata.files(dir(path = output.dir, pattern = "genotype.probs.txt",
                     full.names = TRUE), cross = sampletype)

  # Create a single founder allele probability file.
  if(sampletype == "DOF1") {

    condense.model.probs(path = output.dir, write = paste(output.dir,
                         "founder.probs.Rdata", sep = "/"), model = "full",
                         cross = sampletype)

  } else {

    condense.model.probs(path = output.dir, write = paste(output.dir,
                         "founder.probs.Rdata", sep = "/"), cross = sampletype)

  } # else
 
  # Create genotype plots if the user requested it.
  if(plot) {

    write.genoprob.plots(path = output.dir, snps = snps)

  } # if(plot)

} # calc.genoprob()


################################################################################
# Helper functions.
read.muga.data = function(array = c("gigamuga", "megamuga", "muga"),
                 method = c("intensity", "allele", "both"),
                 sampletype = c("DO", "CC", "HS")) {

  retval = as.list(rep(NA, 6))
  names(retval)= c("snps", "sex", "code", "x", "y", "geno")

  array  = match.arg(array)
  method = match.arg(method)
  sampletype = match.arg(sampletype )

  snpfile  = NULL
  sexfile  = NULL
  codefile = NULL
  xfile    = NULL
  yfile    = NULL
  genofile = NULL

  if(array == "muga") {

    snpfile  = "ftp://ftp.jax.org/MUGA/muga_snps.Rdata"
    sexfile  = "ftp://ftp.jax.org/MUGA/muga_sex.Rdata"
    codefile = "ftp://ftp.jax.org/MUGA/muga_code.Rdata"
    xfile    = "ftp://ftp.jax.org/MUGA/muga_x.Rdata"
    yfile    = "ftp://ftp.jax.org/MUGA/muga_y.Rdata"
    genofile = "ftp://ftp.jax.org/MUGA/muga_geno.Rdata"

    if(sampletype == "HS") { 
      codefile = "ftp://ftp.jax.org/MUGA/muga_hs_code.Rdata"
    } # if(sampletype == "HS")

  } else if(array == "megamuga") {

    snpfile  = "ftp://ftp.jax.org/MUGA/MM_snps.Rdata"
    sexfile  = "ftp://ftp.jax.org/MUGA/MM_sex.Rdata"
    codefile = "ftp://ftp.jax.org/MUGA/MM_code.Rdata"
    xfile    = "ftp://ftp.jax.org/MUGA/MM_x.Rdata"
    yfile    = "ftp://ftp.jax.org/MUGA/MM_y.Rdata"
    genofile = "ftp://ftp.jax.org/MUGA/MM_geno.Rdata"

    if(sampletype == "HS") { 
      codefile = "ftp://ftp.jax.org/MUGA/MM_hs_code.Rdata"
    } # if(sampletype == "HS")

  } else if(array == "gigamuga") {

    snpfile  = "ftp://ftp.jax.org/MUGA/GM_snps.Rdata"
    sexfile  = "ftp://ftp.jax.org/MUGA/GM_sex.Rdata"
    codefile = "ftp://ftp.jax.org/MUGA/GM_code.Rdata"
    xfile    = "ftp://ftp.jax.org/MUGA/GM_x.Rdata"
    yfile    = "ftp://ftp.jax.org/MUGA/GM_y.Rdata"
    genofile = "ftp://ftp.jax.org/MUGA/GM_geno.Rdata"

    if(sampletype == "HS") { 
      codefile = "ftp://ftp.jax.org/MUGA/GM_hs_code.Rdata"
    } # if(sampletype == "HS")

  } else {
    stop(paste("read.muga.data: Unknown array:", array))
  } # else

  if(array == "muga") {

    # Get the MUGA SNPs and subset them to keep the ones that we use for 
    # genotyping.
    # We have to put this line in to satisfy R CMD build --as-cran
    muga_snps = NULL
    load(url(snpfile))
    snps = muga_snps
    snps = snps[snps[,1] %in% muga.snps.to.keep,]
    snps = snps[!is.na(snps[,4]),]

    # We have to put this line in to satisfy R CMD build --as-cran
    muga_sex  = NULL
    muga_code = NULL
    load(url(sexfile))
    load(url(codefile))

    retval[[1]] = snps
    retval[[2]] = muga_sex
    retval[[3]] = muga_code

    if(method == "intensity" | method == "both") {

      muga_x = NULL
      muga_y = NULL
      load(url(xfile))
      load(url(yfile))

      retval[[4]] = t(muga_x)
      retval[[5]] = t(muga_y)

    } # if(method == "intensity" | method == "both")

    if(method == "allele" | method == "both") {

      muga_geno = NULL
      load(url(genofile))

      retval[[6]] = t(muga_geno)

    } # if(method == "allele" | method == "both")

  } else if(array == "megamuga") {

    # Get the MegaMUGA SNPs and subset them to keep the ones that we use for 
    # genotyping.
    # We have to put this line in to satisfy R CMD build --as-cran
    MM_snps = NULL
    load(url(snpfile))
    snps = MM_snps
    snps = snps[!is.na(snps[,4]),]

    # We have to put this line in to satisfy R CMD build --as-cran
    MM_sex  = NULL
    MM_code = NULL
    load(url(sexfile))
    load(url(codefile))

    retval[[1]] = snps
    retval[[2]] = MM_sex
    retval[[3]] = MM_code

    if(method == "intensity" | method == "both") {

      MM_x = NULL
      MM_y = NULL
      load(url(xfile))
      load(url(yfile))

      retval[[4]] = t(MM_x)
      retval[[5]] = t(MM_y)

    } # if(method == "intensities" | method == "both")

    if(method == "allele" | method == "both") {

      muga_geno = NULL
      load(url(genofile))
      retval[[6]] = t(MM_geno)

    } # if(method == "allele" | method == "both")

  } else if(array == "gigamuga") {

    # Get the GigaMUGA SNPs and subset them to keep the ones that we use for 
    # genotyping.
    # We have to put this line in to satisfy R CMD build --as-cran
    GM_snps = NULL
    load(url(snpfile))
    snps = GM_snps
    snps = snps[!is.na(snps[,4]),]
    snps = snps[snps[,3] > 0,]

    # We have to put this line in to satisfy R CMD build --as-cran
    GM_sex  = NULL
    GM_code = NULL
    load(url(sexfile))
    load(url(codefile))

    retval[[1]] = snps
    retval[[2]] = GM_sex
    retval[[3]] = GM_code

    if(method == "intensity" | method == "both") {

      GM_x = NULL
      GM_y = NULL
      load(url(xfile))
      load(url(yfile))

      retval[[4]] = t(GM_x)
      retval[[5]] = t(GM_y)

    } # if(method == "intensities" | method == "both")

    if(method == "allele" | method == "both") {

      GM_geno = NULL
      load(url(genofile))

      retval[[6]] = t(GM_geno)

    } # if(method == "allele" | method == "both")

  } else {
    stop(paste("read.muga.data: Unknown array:", array))
  } # else

  retval[!is.na(retval)]          

} # read.muga.data()


# Synch up the SNPs in snps and the data list passed in.
# snps: data.frame containing 4 columns: marker ID, chr, Mb position, cM position.
# data: The data list that contains either x and y or geno, sex and generation.
# founders: The founders list that contains either x and y or geno, sex and code.
synchronize.snps = function(markers, data, founders) {

  # Figure out which list elements are matrices in data and founders.
  data.wh     = which(!sapply(sapply(data, dim), is.null))
  founders.wh = which(!sapply(sapply(founders, dim), is.null))

  if(length(data.wh) == 0) {
    stop(paste("The data list does not contain any matrices. It must",
         "contain either x and y matrices or a geno matrix."))
  } # if(length(data.wh) == 0)

  if(length(founders.wh) == 0) {
    stop(paste("The founder list does not contain any matrices. It must",
         "contain either x and y matrices or a geno matrix."))
  } # if(length(data.wh) == 0)

  # Get common markers.
  id = markers[,1]
  for(i in 1:length(data.wh)) {
    id = intersect(id, colnames(data[[data.wh[i]]]))
  } # for(i)

  for(i in 1:length(founders.wh)) {
    id = intersect(id, colnames(founders[[founders.wh[i]]]))
  } # for(i)

  markers = markers[markers[,1] %in% id,]
  for(i in 1:length(data.wh)) {
    data[[data.wh[i]]] = data[[data.wh[i]]][,id]
  } # for(i)

  for(i in 1:length(founders.wh)) {
    founders[[founders.wh[i]]] = founders[[founders.wh[i]]][,id]
  } # for(i)

  return(list(markers = markers, data = data, founders = founders))

} # synchronize.snps()


# Keep only the DO founders with codes that are no NA.
# Founders is expected to contain matrices with genotypes or intensties
# that contain samples in rows and markers in columns.
# Founders will also contain vecters 'sex' and 'code'.
keep.do.founders = function(founders) {

  if(all(names(founders) != "code")) {
    stop(paste("keep.do.founders: \'code\' was not found in \'founders\'."))
  } # if(all(names(founders) != "code"))

  vec = which(sapply(founders, length) == length(founders$code))
  mat = sapply(founders, ncol)
  mat = which(!sapply(mat, is.null))

  keep = which(!is.na(founders$code))
  for(i in vec) {
    founders[[i]] = founders[[i]][keep]
  } # for(i)

  for(i in mat) {
    founders[[i]] = founders[[i]][keep,,drop = FALSE]
  } # for(i)

  return(founders)

} # keep.do.founders()

# Verify that the alleles are the same at each marker in the
# founders and data genotypes.
check.alleles = function(data, founders) {

  if(!"geno" %in% names(data)) {
    stop("check.alleles: \'data\' argument does not contain a \'geno\' element.")
  } # if(!"g" %in% names(data))

  if(!"geno" %in% names(founders)) {
    stop("check.alleles: \'founders\' argument does not contain a \'geno\' element.")
  } # if(!"g" %in% names(founders))

  stopifnot(ncol(data$geno) == ncol(founders$geno))

  dg = apply(data$geno, 2, unique)
  fg = apply(founders$geno, 2, unique)

  dg = lapply(dg, sort)
  fg = lapply(fg, sort)

  dg = lapply(dg, function(z) { z[z != "N"] })
  fg = lapply(fg, function(z) { z[z != "N"] })

  if(!all(sort(unique(unlist(dg))) == sort(unique(unlist(fg))))) {
    stop(paste("All alleles do not match in founders and data genotypes.",
         "Please verify that founders$geno and data$geno have the same allele calls."))
  } # if(!all(sort(unique(unlist(dg))) == sort(unique(unlist(fg)))))

} # check.alleles()
