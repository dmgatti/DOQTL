################################################################################
# Check that the dimensions of the data arrays and vectors are consistent.
# Daniel Gatti
# Dan.Gatti@jax.org
# June 20, 2013
################################################################################
check.data.dims = function(data) {

  dims = sapply(data, dim)

  arrays  = which(!is.null(dims))
  vectors = which(is.null(dims))

  for(i in arrays) {
    for(j in vectors) {
      if(length(data[[j]]) != nrow(data[[i]])) {
        stop(paste("Data dimensions not consistent. The number of rows in the",
             names(data)[i], "matrix (", nrow(data[[i]]),
             ")does not equal the length of the", names(data)[j], "vector (",
             length(data[[j]]), ")."))
      } # if(length(data[[j]]) != nrow(data[[i]])
    } # for(j)
  } # for(i)

} # check.data.dims()
