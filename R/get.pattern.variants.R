# Given a VCF class (from VariantAnnotation) and a set of strain names that are
# in the sample names of the VCF, return variants for which the strain.subset
# has the alternate (or different) alleles.
# variants: VCF object of the kind returned by readVcf().
# strain.subset = character vector of strains names that match some strains in
#                 the samples(variants).
get.pattern.variants = function(variants, strain.subset = NULL) {

  if(!is.null(strain.subset)) {
 
    type = attr(variants, "type")

    # Verify that all of the strains in the subset are in the colnames of variants.
    if(!all(strain.subset %in% samples(variants))) {
      stop("all of the subset strains are not in colnames(variants).")
    } # if(!all(strain.subset %in% colnames(variants)))

    # Separate the header from the variants.
    hdr = variants[,1:5]
    variants = variants[,-1:-5]

    # If there are confidence values, strip them out.
    variants = strip.quality.columns(variants)

    # Make subset indices.
    ss = which(colnames(variants)  %in% strain.subset)
    mm = which(!colnames(variants) %in% strain.subset)

    # First, get all of the variants where the members of each subset have the
    # same allele.
    ss.allele = apply(variants[,ss, drop = FALSE], 1, unique)
    mm.allele = apply(variants[,mm, drop = FALSE], 1, unique)
    ss.allele = sapply(ss.allele, function(a) { a[!is.na(a)] })
    mm.allele = sapply(mm.allele, function(a) { a[!is.na(a)] })

    # Keep only the variants where each subset has only one allele.
    keep = which((sapply(ss.allele, length) == 1) &
                 (sapply(mm.allele, length) == 1))
    hdr  = hdr[keep,,drop = FALSE]
    variants = variants[keep,,drop = FALSE]
    ss.allele = unlist(ss.allele[keep])
    mm.allele = unlist(mm.allele[keep])

    # Keep only the variants that differ between groups.
    keep = ss.allele != mm.allele
    hdr  = hdr[keep,,drop = FALSE]
    variants = variants[keep,,drop = FALSE]
    variants = data.frame(hdr, variants, stringsAsFactors = FALSE)
    attr(variants, "type") = type

  } # if(!is.null(strain.subset))

  return(variants)

} # get.pattern.variants()

