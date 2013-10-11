################################################################################
# Convert the genotype data from A,C,G,T format to A, H, B, N.
# Daniel Gatti
# Dan.Gatti@jax.org
# Arguments: geno: character matrix with genotypes. Sample IDs in rownames
#                  & SNP IDs in colnames. num.samples x num.SNPs.
################################################################################
convert.genotypes = function(geno) {

  # Get the alleles for each SNP.
  tbl = apply(geno, 2, table)
  # Remove the no-calls for the moment.  We'll deal with them below.
  tbl = lapply(tbl, function(a) { a[names(a) != "-"] })
  # Count the number of alleles.
  tbl.len = sapply(tbl, length)
  # If we do not have more than 3 alleles, stop.
  if(max(tbl.len) > 4) {
    stop("The number of alleles at some SNPs is greater than 4.")
  } # if(max(tbl.len) > 4)

  # Process the single allele genotypes first.
  ss = which(tbl.len == 1)
  geno[,ss] = matrix("A", nrow(geno), length(ss), byrow = T)

  # Process two allele genotypes.
  ss = which(tbl.len == 2)
  alleles = sapply(tbl[ss], names)

  for(i in 1:ncol(alleles)) {
    if(alleles[1,i] %in% c("AA", "CC", "GG", "TT")) {
      geno[geno[,ss[i]] == alleles[1,i], ss[i]] = "A"
      if(alleles[2,i] %in% c("AA", "CC", "GG", "TT")) {
        geno[geno[,ss[i]] == alleles[2,i], ss[i]] = "B"
      } else {
        geno[geno[,ss[i]] == alleles[2,i], ss[i]] = "H"
      } # else
    } else {
      geno[geno[,ss[i]] == alleles[1,i], ss[i]] = "H"
      geno[geno[,ss[i]] == alleles[2,i], ss[i]] = "A"
    } # else
  } # for(i)

  # Process three allele genotypes.
  ss = which(tbl.len == 3)
  alleles = sapply(tbl[ss], names)

  for(i in 1:ncol(alleles)) {
    geno[geno[,ss[i]] == alleles[1,i], ss[i]] = "A"
    geno[geno[,ss[i]] == alleles[2,i], ss[i]] = "H"
    geno[geno[,ss[i]] == alleles[3,i], ss[i]] = "B"
  } # for(i)

  # Replace the missing data symbol with = "N".
  geno[geno == "-"] = "N"

  return(geno)

} # convert.genotypes()