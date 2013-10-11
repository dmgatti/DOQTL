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
#                   samples. The muga and megamuga opoitns are for the Mouse
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
                         plot = T, array = c("megamuga", "muga", "other"),
                         sampletype = c("DO", "CC", "DOF1", "other"), method = c("intensity",
                         "allele"), founders, transprobs, snps) {

  if(missing(data)) {
    stop(paste("data is missing. Please supply data to process."))
  } # if(missing(data))

  method = match.arg(method)
  sampletype = match.arg(sampletype)
  array = match.arg(array)
  attr(data, "method") = method
  attr(data, "array")  = array
  attr(data, "sampletype") = sampletype

  # If the sampletype is CC or DO and array is muga or megamuga, then get the 
  # founder data.
  if(sampletype == "CC" | sampletype == "DO") {
    if(array == "muga") {

      # Get the MUGA SNPs and subset them to keep the ones that we use for 
      # genotyping.
      load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
      snps = muga_snps

      data(muga.snps.to.keep)
      snps = snps[snps[,1] %in% muga.snps.to.keep,]
      snps = snps[!is.na(snps[,4]),]

      load(url("ftp://ftp.jax.org/MUGA/muga_sex.Rdata"))
      load(url("ftp://ftp.jax.org/MUGA/muga_code.Rdata"))

      if(method == "allele") {
        load(url("ftp://ftp.jax.org/MUGA/muga_geno.Rdata"))
        founders = list(geno = t(muga_geno[rownames(muga_geno) %in% snps[,1],]),
                        sex = muga_sex, code = muga_code,
                        states = create.genotype.states(LETTERS[1:8]))
        if(!"geno" %in% names(data)) {
          stop(paste("There is no item called 'geno' in the data list. Please",
               "add the allele calls in a matrix called 'geno' when method",
               "= 'alelle'."))
        } else {
          data$geno = as.matrix(data$geno[,colnames(data$geno) %in% snps[,1]])
          snps = snps[snps[,1] %in% colnames(data$geno),]
        } # else

      } else {
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
          data$x = as.matrix(data$x[,colnames(data$x) %in% snps[,1]])
          data$y = as.matrix(data$y[,colnames(data$y) %in% snps[,1]])
          snps = snps[snps[,1] %in% colnames(data$x),]
        } # else

      } # else

    } else if(array == "megamuga") {
      
      # Get the SNPs on the MegaMUGA and subset them to keep the ones for
      # the CC/DO mice.
      load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
      snps = MM_snps

      # Keep only the collaborative cross SNPs.
#      snps = snps[snps$Collaborative.Cross == 1,]
      snps = snps[!is.na(snps[,4]),]
      snps = snps[,1:4]

      load(url("ftp://ftp.jax.org/MUGA/MM_sex.Rdata"))
      load(url("ftp://ftp.jax.org/MUGA/MM_code.Rdata"))

      if(method == "allele") {

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

        data$geno = as.matrix(data$geno[,colnames(data$geno) %in% snps[,1]])
        snps = snps[snps[,1] %in% colnames(data$geno),]
        rm(MM_geno)

      } else if(method == "intensity") {

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
        data$x = as.matrix(data$x[,colnames(data$x) %in% snps[,1]])
        data$y = as.matrix(data$y[,colnames(data$y) %in% snps[,1]])
        snps = snps[snps[,1] %in% colnames(data$x),]
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
  } else if (smapletype == "DOF1") {
    stop("DOF1 samples not implemented yet.")
    # In this case, the user will need to supply, at a minimum, the 9 founders.
    # We may be able to supply the 8 DO founders and only require the 9th founder.
  } else {
    stop("Non-CC and non-DO samples not implemented yet.")
    # The user is going to have to supply everything: SNPs, founders & trans.probs.
  } # else

  # Set sex identifiers to upper case.
  data$sex = as.matrix(data$sex)
  data$sex = toupper(data$sex)

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
    trans.prob.fxn = generic.trans.probs
  } else if(sampletype == "other") {
    trans.prob.fxn = generic.trans.probs
  } # else

  chr = as.character(chr)

  # Verify that SNPs are sorted.
  tmp = split(snps[,1:4], snps[,2])
  tmp = lapply(tmp, function(z) { diff(z[,4]) })
  if(any(unlist(sapply(tmp, "<", 0)), na.rm = T)) {
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
     calc.genoprob.alleles(data = data, chr = chr, founders = founders, snps = snps,
                           output.dir = output.dir, trans.prob.fxn = trans.prob.fxn,
                           plot = plot)

  } else if(method == "intensity") {

     # Convert the sample names into names that can be written out 
     # as files.
     rownames(data$x) = make.names(rownames(data$x))
     rownames(data$y) = make.names(rownames(data$y))
     names(data$sex) = make.names(names(data$sex))
     if(attr(data, "sampletype") == "DO") {
       names(data$gen) = make.names(names(data$gen))
     } # if(attr(data, "sampletype") == "DO")

     calc.genoprob.intensity(data = data, chr = chr, founders = founders, 
                             snps = snps, output.dir = output.dir, 
                             trans.prob.fxn = trans.prob.fxn, plot = plot)
  } # else if(method = "intensity")

  # Convert the *.genotype.probs.txt files to *.Rdata files.
  create.Rdata.files(dir(path = output.dir, pattern = "genotype.probs.txt",
                     full.names = T))

  # Create a single founder allele probability file.
  condense.model.probs(path = output.dir, write = paste(output.dir,
                       "founder.probs.Rdata", sep = "/"))
 
  # Create genotype plots if the user requested it.
  if(plot) {
#    write.genotype.plots(path = output.dir, snps)
  } # if(plot)
} # calc.genoprob()

