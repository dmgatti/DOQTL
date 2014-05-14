################################################################################
# Given the founder data, a chromosome and the SNPs, produce emission 
# probabilities for the given chromosome.
# Arguments: founders: list containing founder genotypes and sex.
#            chr: character (or numeric) with a single chromosome.
#            snps: data.frame with 4 columns containing SNPs. SNP IDs, Chr,
#                  Mb & cM in columns 1:4.
#            sex: chracter, either "F" or "M". Only used on sex chromosomes.
# Daniel Gatti
# Dan.Gatti@jax.org
# April, 17, 2013
################################################################################
emission.probs.allele = function(founders, chr, snps, sex = c("F", "M")) {
  print("Getting emission probabilities from founder data...")
  if(ncol(founders$geno) != nrow(snps)) {
    stop(paste("emission.probs.alleles: The number of columns in founders$geno",
         "is not equal to the number of SNPs."))
  } # if(ncol(founders$geno) != nrow(snps))
  if(any(snps[,1] != colnames(founders$geno))) {
    stop(paste("emission.probs.alleles: The SNP IDs in snps do not match those",
         "in geno. Please verify that the SNP IDs are identical."))
  } # if(any(snps[,1] != colnames(founders$geno)))
  # Convert genotypes to numeric values (0 = A, 1 = H, 2 = B, 3 = N)
  chr = as.character(chr[1])
  geno = convert.allele.calls(as.matrix(founders$geno))
  symbols = sort(unique(as.vector(geno)))
  retval = NULL
  # Autosomes.
  if(!is.na(as.numeric(chr))) {
    print(paste("Found", nrow(snps), "SNPs on Chr", chr))
    print(paste("Found", length(founders$states), "genotype states."))
    retval = tabulate.geno(geno = geno, founders = founders,
             symbols = symbols)
  } else if(chr == "X") {
    sex = match.arg(sex)
    sample.subset = which(founders$sex == sex)
    print(paste("Found", nrow(snps), "SNPs on Chr", chr))
    print(paste("Found", length(founders$states), "genotype states."))
    if(sex == "F") {
      retval = tabulate.geno(geno = geno[sample.subset,], founders = founders,
               symbols = symbols)
    } else if (sex == "M") {
      founders$states = paste(founders$states, founders$states, sep = "")
      retval = tabulate.geno(geno = geno[sample.subset,], founders = founders,
               symbols = symbols)
    } else {
      stop("emission.probs.allele: Unknown sex.")
    } # else
  } else if(chr == "Y") {
    print("Male Y Chr...")
    male = which(founders$sex == "M")
    print(paste("Found", length(founders$states), "genotype states."))
    founders$states = paste(founders$states, founders$states, sep = "")
    retval = tabulate.geno(geno = geno[male,], founders = founders, 
             symbols = symbols)
  } else if(chr == "M") {
    print("M Chr...")
    print(paste("Found", length(founders$states), "genotype states."))
    founders$states = paste(founders$states, founders$states, sep = "")
    retval = tabulate.geno(geno = geno, founders = founders, 
             symbols = symbols)
  } # else if(chr == "M")
  return(log(retval))
} # emission.probs.allele()

# Helper function to assign the genotype probabilities.
tabulate.geno = function(geno, founders, symbols) {
    retval = array(rnorm(n = length(symbols) * length(founders$states) * ncol(geno),
          mean = 0.01, sd = 0.0001),
	      c(length(symbols), length(founders$states), ncol(geno)),
          dimnames = list(symbols, founders$states, colnames(geno)))
    for(s in founders$states) {
      # Get the rows where this genotype occurs.
      rows = which(founders$code == s)
      curr.geno = data.frame(geno[rows,, drop = FALSE], row.names = NULL)
      # Convert each column (SNP) into a factor with four levels, one for 
      # each symbol.
      curr.geno = lapply(curr.geno, factor, levels = symbols)
      # Tabulate the occurrence of each symbol.
      tbl = lapply(curr.geno, table)
      tbl = lapply(tbl, function(a) { a / length(rows) })
      tmp = matrix(unlist(tbl), dim(retval)[1], dim(retval)[3])
      eq.zero = which(tmp == 0)
      tmp[eq.zero] = rnorm(n = length(eq.zero), mean = 0.01, sd = 0.0001)
      retval[,s,] = tmp / matrix(colSums(tmp),
                    nrow = nrow(tmp), ncol = ncol(tmp), byrow = TRUE)
    } # for(s)	
    return(retval)
} # tabulate.geno()
