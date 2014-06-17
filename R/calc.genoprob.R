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
#                        'CC', 'DOF1' or 'other'.
#            method: character, that inidicates the method of genotype reconstruction.
#                    'intensity' uses X and Y array intensities and 'allele' uses
#                    the genotype calls.
#            founders: list, containing founder data. Must contain allele calls
#                         or X & Y intensites, founder names, codes and sexes.
#            transprobs: list, containing transition probabilities between each pair
#                        of SNPs on each chromosome.
#            snps: data.frame with SNP IDs, chr, Mb and cM positions. For custom
#                  arrays.
calc.genoprob = function(data, chr = "all", output.dir = ".",
                         plot = TRUE, array = c("megamuga", "muga", "other"),
                         sampletype = c("DO", "CC", "DOF1", "other"), method =
                         c("intensity", "allele"), founders, transprobs, snps) {

  if(missing(data)) {
    stop(paste("data is missing. Please supply data to process."))
  } # if(missing(data))

  if(!is.null(data$sex)) {
    if(is.matrix(data$sex) || is.list(data$sex)) { 
      nm = rownames(data$sex)
      data$sex = as.character(data$sex[,1])
      names(data$sex) = nm
    }
  } else {
    stop("data$sex cannot be null, even if all animals are the same sex.")
  } # else

  if(!is.null(data$gen)) {
    if(is.matrix(data$gen) || is.list(data$gen)) { 
      nm = rownames(data$gen)
      data$gen = as.character(data$gen[,1])
      names(data$gen) = nm
    } # if(is.matrix(data$gen) || is.list(data$gen))
  } # if(!is.null(data$gen))

  method = match.arg(method)
  sampletype = match.arg(sampletype)
  array = match.arg(array)
  attr(data, "method") = method
  attr(data, "array")  = array
  attr(data, "sampletype") = sampletype

  # If the sampletype is CC or DO and array is muga or megamuga, then get the 
  # founder data.
  ### CC or DO ###
  if(sampletype == "CC" | sampletype == "DO") {
    ### MUGA ###
    if(array == "muga") {
      # Get the MUGA SNPs and subset them to keep the ones that we use for 
      # genotyping.
	# We have to put this line in to satisfy R CMD build --as-cran
      muga_snps = NULL
      load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
      snps = muga_snps
      snps = snps[snps[,1] %in% muga.snps.to.keep,]
      snps = snps[!is.na(snps[,4]),]
      # We have to put this line in to satisfy R CMD build --as-cran
      muga_sex  = NULL
      muga_code = NULL
      load(url("ftp://ftp.jax.org/MUGA/muga_sex.Rdata"))
      load(url("ftp://ftp.jax.org/MUGA/muga_code.Rdata"))
      
      ### Allele Call ###
      if(method == "allele") {
        # We have to put this line in to satisfy R CMD build --as-cran
        muga_geno = NULL
        load(url("ftp://ftp.jax.org/MUGA/muga_geno.Rdata"))
        founders = list(geno = t(muga_geno[rownames(muga_geno) %in% snps[,1],]),
                        sex = muga_sex, code = muga_code,
                        states = create.genotype.states(LETTERS[1:8]))
        if(!"geno" %in% names(data)) {
          stop(paste("There is no item called 'geno' in the data list. Please",
               "add the allele calls in a matrix called 'geno' when method",
               "= 'alelle'."))
        } else {
          # Align the SNPs in the data and founders.
          keep = intersect(intersect(snps[,1], colnames(data$geno)), colnames(founders$geno))
          data$geno = as.matrix(data$geno[,colnames(data$geno) %in% keep])
          founders$geno = as.matrix(founders$geno[,colnames(founders$geno) %in% keep])
          snps = snps[snps[,1] %in% keep,]
          data$geno = data$geno[,match(snps[,1], colnames(data$geno))]
          founders$geno = founders$geno[,match(snps[,1], colnames(founders$geno))]
        } # else
     ### Intensity ###
      } else {

        # We have to put this line in to satisfy R CMD build --as-cran
        muga_x = NULL
        muga_y = NULL

        load(url("ftp://ftp.jax.org/MUGA/muga_x.Rdata"))
        load(url("ftp://ftp.jax.org/MUGA/muga_y.Rdata"))
        founders = list(x = t(muga_x[rownames(muga_x) %in% snps[,1],]),
                        y = t(muga_y[rownames(muga_y) %in% snps[,1],]), 
                        sex = muga_sex, code = muga_code,
                        states = create.genotype.states(LETTERS[1:8]))

        if(!"x" %in% names(data) | !"y" %in% names(data)) {
          stop(paste("There 'x' or 'y' matrices are missing from the data list.",
               "Please add the X and Y intensities in matrices called 'x' and 'y'",
               "when method = 'intensity'."))
        } else {
          # Align the SNPs in the data and founders.
          keep = intersect(intersect(snps[,1], colnames(data$x)), colnames(founders$x))
          data$x = as.matrix(data$x[,colnames(data$x) %in% keep])
          data$y = as.matrix(data$y[,colnames(data$y) %in% keep])
          founders$x = as.matrix(founders$x[,colnames(founders$x) %in% keep])
          founders$y = as.matrix(founders$y[,colnames(founders$y) %in% keep])
          snps = snps[snps[,1] %in% keep,]
          data$x = data$x[,match(snps[,1], colnames(data$x))]
          data$y = data$y[,match(snps[,1], colnames(data$y))]
          founders$x = founders$x[,match(snps[,1], colnames(founders$x))]
          founders$y = founders$y[,match(snps[,1], colnames(founders$y))]
        } # else
      } # else

    ### MegaMUGA ###
    } else if(array == "megamuga") {
      
      # Get the SNPs on the MegaMUGA and subset them to keep the ones for
      # the CC/DO mice.
      # We have to put this line in to satisfy R CMD build --as-cran
      MM_snps = NULL
      load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
      snps = MM_snps

      # Keep only the collaborative cross SNPs.
      snps = snps[snps$Collaborative.Cross == 1 | snps$MUGA == 1 |
	              snps$C57BL.6 == 1,]
      snps = snps[!is.na(snps[,4]),]
      snps = snps[,1:4]

      # We have to put this line in to satisfy R CMD build --as-cran
      MM_sex = NULL
      MM_code = NULL
      load(url("ftp://ftp.jax.org/MUGA/MM_sex.Rdata"))
      load(url("ftp://ftp.jax.org/MUGA/MM_code.Rdata"))
 
     ### Allelle Call ###
      if(method == "allele") {
        # We have to put this line in to satisfy R CMD build --as-cran
        MM_geno = NULL
        load(url("ftp://ftp.jax.org/MUGA/MM_geno.Rdata"))
        # Remove founders with too much missing data.
        remove = which(colSums(MM_geno == "N") > 5000)
        MM_geno = MM_geno[,-remove]
        MM_sex  = MM_sex[-remove]
        MM_code = MM_code[-remove]
        founders = list(geno = t(MM_geno[rownames(MM_geno) %in% snps[,1],]),
                        sex = MM_sex, code = MM_code,
                        states = create.genotype.states(LETTERS[1:8]))
        keep = which(!is.na(MM_code))
        founders$geno = founders$geno[keep,]
        founders$sex  = founders$sex[keep]
        founders$code = founders$code[keep]
        # Align the SNPs in the data and founders.
        keep = intersect(intersect(snps[,1], colnames(data$geno)), colnames(founders$geno))
        data$geno = as.matrix(data$geno[,colnames(data$geno) %in% keep])
        founders$geno = as.matrix(founders$geno[,colnames(founders$geno) %in% keep])
        snps = snps[snps[,1] %in% keep,]
        data$geno = data$geno[,match(snps[,1], colnames(data$geno))]
        founders$geno = founders$geno[,match(snps[,1], colnames(founders$geno))]
        data$geno = as.matrix(data$geno[,colnames(data$geno) %in% snps[,1]])
        snps = snps[snps[,1] %in% colnames(data$geno),]
        rm(MM_geno)
      ### Intensity ###
      } else if(method == "intensity") {
        # We have to put this line in to satisfy R CMD build --as-cran
        MM_x = NULL
        MM_y = NULL
        load(url("ftp://ftp.jax.org/MUGA/MM_x.Rdata"))
        load(url("ftp://ftp.jax.org/MUGA/MM_y.Rdata"))
        # Remove founders with too much missing data.
        remove = which(colSums(MM_x < 0) > 1000)
        MM_x = MM_x[,-remove]
        MM_y = MM_y[,-remove]
        MM_sex  = MM_sex[-remove]
        MM_code = MM_code[-remove]
        founders = list(x = t(MM_x[rownames(MM_x) %in% snps[,1],]), 
                        y = t(MM_y[rownames(MM_y) %in% snps[,1],]),
                        sex = MM_sex, code = MM_code,
                        states = create.genotype.states(LETTERS[1:8]))
        keep = which(!is.na(MM_code))
        founders$x = founders$x[keep,]
        founders$y = founders$y[keep,]
        founders$sex  = founders$sex[keep]
        founders$code = founders$code[keep]
        # Align the SNPs in the data and founders.
        keep = intersect(intersect(snps[,1], colnames(data$x)), colnames(founders$x))
        data$x = as.matrix(data$x[,colnames(data$x) %in% keep])
        data$y = as.matrix(data$y[,colnames(data$y) %in% keep])
        founders$x = as.matrix(founders$x[,colnames(founders$x) %in% keep])
        founders$y = as.matrix(founders$y[,colnames(founders$y) %in% keep])
        snps = snps[snps[,1] %in% keep,]
        data$x = data$x[,match(snps[,1], colnames(data$x))]
        data$y = data$y[,match(snps[,1], colnames(data$y))]
        founders$x = founders$x[,match(snps[,1], colnames(founders$x))]
        founders$y = founders$y[,match(snps[,1], colnames(founders$y))]
        if(!all(colnames(data$x) == colnames(founders$x))) {
          stop("SNP names in data$x and founders$x do not match for MegaMUGA.")
        } # if(!all(colnames(data$x) == colnames(founder$x)))
        if(!all(colnames(data$y) == colnames(founders$y))) {
          stop("SNP names in data$y and founders$y do not match for MegaMUGA.")
        } # if(!all(colnames(data$y) == colnames(founder$y)))
        
        rm(MM_x, MM_y)
      } # else
      rm(MM_sex, MM_code)
      attr(founders, "method") = method

    ### Custome Array ###
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
        # Synch up the SNPs in the data.
        founders$geno = founders$geno[,colnames(founders$geno) %in% snps[,1]]
        data$geno = data$geno[,colnames(data$geno) %in% snps[,1]]
        snps = snps[snps[,1] %in% colnames(founders$geno),]
        founders$geno = as.matrix(founders$geno)
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
        # Synch up the SNPs in the data.
        founders$x = as.matrix(founders$x[,colnames(founders$x) %in% snps[,1]])
        founders$y = as.matrix(founders$y[,colnames(founders$y) %in% snps[,1]])
        data$x = as.matrix(data$x[,colnames(data$x) %in% snps[,1]])
        data$y = as.matrix(data$y[,colnames(data$y) %in% snps[,1]])
        snps = snps[snps[,1] %in% colnames(founders$x),]
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
  ### DO/F1 ###
  } else if (sampletype == "DOF1") {
    ### MegaMUGA ###
    if(array == "megamuga" ) {
      # Get the SNPs on the MegaMUGA and subset them to keep the ones for
      # the CC/DO mice.
      # We have to put this line in to satisfy R CMD build --as-cran
      MM_snps = NULL
      load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
      snps = MM_snps
      # Keep only the collaborative cross SNPs.
      snps = snps[snps$Collaborative.Cross == 1 | snps$MUGA == 1 |
	              snps$C57BL.6 == 1,]
      snps = snps[!is.na(snps[,4]),]
      snps = snps[,1:4]
      # We have to put this line in to satisfy R CMD build --as-cran
      MM_sex = NULL
      MM_code = NULL
      load(url("ftp://ftp.jax.org/MUGA/MM_sex.Rdata"))
      load(url("ftp://ftp.jax.org/MUGA/MM_code.Rdata"))
      ### Allele Call ###
      if(method == "allele") {
        # We have to put this line in to satisfy R CMD build --as-cran
        MM_geno = NULL
        load(url("ftp://ftp.jax.org/MUGA/MM_geno.Rdata"))
        # Remove founders with too much missing data.
        remove = which(colSums(MM_geno == "N") > 5000)
        MM_geno = MM_geno[,-remove]
        MM_sex  = MM_sex[-remove]
        MM_code = MM_code[-remove]
        founders = list(geno = t(MM_geno[rownames(MM_geno) %in% snps[,1],]),
                        sex = MM_sex, code = MM_code,
                        states = create.genotype.states(LETTERS[1:8]),
                        direction = founders$direction)
        keep = which(!is.na(MM_code))
        founders$geno = founders$geno[keep,]
        founders$sex  = founders$sex[keep]
        founders$code = founders$code[keep]
        # Align the SNPs in the data and founders.
        keep = intersect(intersect(snps[,1], colnames(data$geno)), colnames(founders$geno))
        data$geno = as.matrix(data$geno[,colnames(data$geno) %in% keep])
        founders$geno = as.matrix(founders$geno[,colnames(founders$geno) %in% keep])
        snps = snps[snps[,1] %in% keep,]
        data$geno = data$geno[,match(snps[,1], colnames(data$geno))]
        founders$geno = founders$geno[,match(snps[,1], colnames(founders$geno))]
        data$geno = as.matrix(data$geno[,colnames(data$geno) %in% snps[,1]])
        snps = snps[snps[,1] %in% colnames(data$geno),]
        rm(MM_geno)
      ### Intensity ###
      } else if(method == "intensity") {
        # We have to put this line in to satisfy R CMD build --as-cran
        MM_x = NULL
        MM_y = NULL
        load(url("ftp://ftp.jax.org/MUGA/MM_x.Rdata"))
        load(url("ftp://ftp.jax.org/MUGA/MM_y.Rdata"))
        # Transpose the MegaMUGA founder data.
        MM_x = t(MM_x)
        MM_y = t(MM_y)
        # Remove founders with too much missing data.
        remove = which(rowSums(MM_x < 0) > 1000)
        MM_x = MM_x[-remove,]
        MM_y = MM_y[-remove,]
        MM_sex  = MM_sex[-remove]
        MM_code = MM_code[-remove]
        # Remove founders with NA founder codes (i.e. non-DO founders)
        remove = which(is.na(MM_code))
        MM_x = MM_x[-remove,]
        MM_y = MM_y[-remove,]
        MM_sex  = MM_sex[-remove]
        MM_code = MM_code[-remove]
        # Keep only the inbred founders and discard the F1s.
        MM_code = MM_code[MM_code %in% paste(LETTERS[1:8], LETTERS[1:8], sep = "")]
        MM_x = MM_x[rownames(MM_x) %in% names(MM_code),]
        MM_y = MM_y[rownames(MM_y) %in% names(MM_code),]
        MM_sex = MM_sex[names(MM_sex) %in% names(MM_code)]
        # Synch up the markers in both founder data sets.
        markers = intersect(colnames(MM_x), colnames(founders$x))
        MM_x = MM_x[,colnames(MM_x) %in% markers]
        MM_y = MM_y[,colnames(MM_y) %in% markers]
        founders$x = founders$x[,colnames(founders$x) %in% markers]
        founders$y = founders$y[,colnames(founders$y) %in% markers]
        rm(markers)
        stopifnot(all(colnames(MM_x) == colnames(founders$x)))
        stopifnot(all(colnames(MM_y) == colnames(founders$y)))
        # Create a new founders list with the 9th founder added.
        x = rbind(MM_x, founders$x)
        y = rbind(MM_y, founders$y)
        sex  = c(MM_sex,  founders$sex)
        code = c(MM_code, founders$code)
        states = paste(LETTERS[1:8], "I", sep = "")
        if(founders$direction == "DOxMUT") {
          states = list(auto = states, X = list(F = states, M = LETTERS[1:8]),
                   founders = LETTERS[1:9])
        } else if(founders$direction == "MUTxDO") {
          states = list(auto = states, X = list(F = states, M = "I"),
                   founders = LETTERS[1:9])
        } else {
          stop(paste("founders$direction =", founders$direction, 
               ". founders$direction must be either MUTxDO or DOxMUT. See",
               "help(calc.genoprob) for more information."))
        } # else
        founders = list(x = x, y = y, sex = sex, code = code,
                        states = states, direction = founders$direction)
        rm(x, y, sex, code, states)
        # Align the SNPs in the data and founders.
        keep = intersect(intersect(snps[,1], colnames(data$x)), colnames(founders$x))
        data$x = as.matrix(data$x[,colnames(data$x) %in% keep])
        data$y = as.matrix(data$y[,colnames(data$y) %in% keep])
        founders$x = as.matrix(founders$x[,colnames(founders$x) %in% keep])
        founders$y = as.matrix(founders$y[,colnames(founders$y) %in% keep])
        snps = snps[snps[,1] %in% keep,]
        data$x = data$x[,match(snps[,1], colnames(data$x))]
        data$y = data$y[,match(snps[,1], colnames(data$y))]
        founders$x = founders$x[,match(snps[,1], colnames(founders$x))]
        founders$y = founders$y[,match(snps[,1], colnames(founders$y))]
        if(!all(colnames(data$x) == colnames(founders$x))) {
          stop("SNP names in data$x and founders$x do not match for MegaMUGA.")
        } # if(!all(colnames(data$x) == colnames(founder$x)))
        if(!all(colnames(data$y) == colnames(founders$y))) {
          stop("SNP names in data$y and founders$y do not match for MegaMUGA.")
        } # if(!all(colnames(data$y) == colnames(founder$y)))
        
        rm(MM_x, MM_y)
      } # else
      rm(MM_sex, MM_code)
      attr(founders, "method") = method
    } else {
      stop(paste("Unsupported array platform:", array))
    } # else
  } else {
    stop("Non-CC and non-DO samples not implemented yet.")
    # The user is going to have to supply everything: SNPs, founders & trans.probs.
  } # else

  # Set sex identifiers to upper case.
  data$sex = as.matrix(data$sex)
  data$sex = toupper(data$sex)
  founders$sex = toupper(founders$sex)

  # Fill in the SNP cM values for which consecutive SNPs have the same cM value,
  # even though they are at different Mb positions.
#  snps = fill.in.snps(snps)

  # Fill in any missing F1s.
  attr(founders, "method") = attr(data, "method")
  founders = add.missing.F1s(founders, snps)

  # Add a slash to the output directory, if required.
  output.dir = add.slash(output.dir)

  # Set the transition probability function.
  trans.prob.fxn = do.trans.probs

  if(sampletype == "CC") {
    trans.prob.fxn = cc.trans.probs
  } else if(sampletype == "DOF1") {
    trans.prob.fxn = dof1.trans.probs
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

  if(method == "allele") {

     # Convert the sample names into names that can be written out 
     # as files.
     rownames(data$geno) = make.names(rownames(data$geno))
     names(data$sex) = make.names(names(data$sex))

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
     snps = split(snps, snps[,2])
     old.warn = options("warn")$warn
     options(warn = -1)
     snps = snps[order(as.numeric(names(snps)))]
     options(warn = old.warn)
     gc()
     calc.genoprob.alleles(data = data, chr = chr, founders = founders, snps = snps,
                           output.dir = output.dir, trans.prob.fxn = trans.prob.fxn,
                           plot = plot)
  } else if(method == "intensity") {

     # Convert the sample names into names that can be written out 
     # as files.
     rownames(data$x) = make.names(rownames(data$x))
     rownames(data$y) = make.names(rownames(data$y))
     names(data$sex) = make.names(names(data$sex))
     if(attr(data, "sampletype") %in% c("DO", "DOF1")) {
       names(data$gen) = make.names(names(data$gen))
     } # if(attr(data, "sampletype") == "DO")

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
     snps = split(snps, snps[,2])
     old.warn = options("warn")$warn
     options(warn = -1)
     snps = snps[order(as.numeric(names(snps)))]
     options(warn = old.warn)
     gc()
     calc.genoprob.intensity(data = data, chr = chr, founders = founders, 
                             snps = snps, output.dir = output.dir, 
                             trans.prob.fxn = trans.prob.fxn, plot = plot)
  } # else if(method = "intensity")
  # Convert the *.genotype.probs.txt files to *.Rdata files.
  create.Rdata.files(dir(path = output.dir, pattern = "genotype.probs.txt",
                     full.names = TRUE))
  # Create a single founder allele probability file.
  if(sampletype == "DOF1") {
    condense.model.probs(path = output.dir, write = paste(output.dir,
                         "founder.probs.Rdata", sep = "/"), model = "full")
  } else {
    condense.model.probs(path = output.dir, write = paste(output.dir,
                         "founder.probs.Rdata", sep = "/"))
  } # else
 
  # Create genotype plots if the user requested it.
  if(plot) {
    load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
    write.genoprob.plots(path = output.dir, snps = muga_snps)
  } # if(plot)
} # calc.genoprob()
