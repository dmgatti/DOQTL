filter.geno.probs <-
function(geno) {

  if(is.null(geno)) {
    stop(paste("filter.geno.probs: geno cannot be null."))
  } # if(is.null(geno))

  if(any(dim(geno) == 0)) {
    stop(paste("filter.geno.probs: one of the dimensions of geno is 0.",
         paste(dim(geno), collapse = ",")))
  } # if(is.null(geno))

  # Take the exp(geno) if the values are not between 0 and 1.
  if(min(geno) < 0) {
    geno = exp(geno)
  } # if(min(geno) < 0)

  if(length(dim(geno)) != 3) {
    stop(paste("filter.geno.probs: geno must be a 3 dimensional array of",
               "genotype probabilities."))
  } # if(is.null(geno))

  # Get the maximum probability for each state at each SNP.
  min.p = apply(geno, 2:3, max)

  # Get the minimum maximum probability for each state across all SNPs.
  lt.45 = min.p < 0.45
  lowp.snps = which(apply(lt.45, 2, any))

  # Remove these SNPs and return.
  if(length(lowp.snps) > 0) {
    geno = geno[,,-lowp.snps]
  } # if(length(lowp.snps) > 0)

  return(geno)
}
