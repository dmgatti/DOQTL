################################################################################
# Convert a genotype matrix from {ACGT} allele calls to 0,1,2,3.
# 0: homozygous A
# 1: heterozygous
# 2: homozygous B
# 3: no call
# Daniel Gatti
# Dan.Gatti@jax.org
# April 8, 2013
################################################################################
convert.allele.calls = function(geno1, geno2) {
  geno = rbind(geno1, geno2)
  # Get the number of alleles at each SNP.
  tbl = apply(geno, 2, unique)
  tbl = lapply(tbl, function(a) { a[a != "H" & a != "N"] })
  tbl = lapply(tbl, sort)
  t2 = unlist(lapply(tbl, paste, collapse = ""))
  t3 = table(t2)
 
  # Create a number genotype matrix.
  g2 = matrix(0, nrow(geno), ncol(geno), dimnames = dimnames(geno))
  # Replace SNPs with a single homozygote with 0.
  t3 = t3[names(t3) != ""]
  single = which(nchar(names(t3)) == 1)
  for(i in single) {
    rng = which(t2 == names(t3)[i])
    g2[,rng][geno[,rng] == names(t3)[i]] = 0
  } # for(i)
  two = which(nchar(names(t3)) == 2)
  for(i in two) {
    rng = which(t2 == names(t3)[i])
    alleles = strsplit(names(t3)[i], split = "")[[1]]
    g2[,rng][geno[,rng] == alleles[1]] = 0
    g2[,rng][geno[,rng] == alleles[2]] = 2
  } # for(i)
  # Replace hets and no calls.
  g2[geno == "H"] = 1
  g2[geno == "N"] = 3
  return(list(g2[1:nrow(geno1),], g2[(nrow(geno1)+1):nrow(g2),]))
} # convert.allele.calls()
