################################################################################
# Impute the genoprobs of a matrix from one set of markers to another.
# Daniel Gatti
# Dan.Gatti@jax.org
# Dec. 20, 2014
################################################################################
# Arguments: data: numeric matrix containing the genotype or haplotype 
#                  probabilities. Named markers in rows and genotypes in 
#                  columns.
#            from: marker information (names, chr, pos) for the markers to 
#                  impute from. Must match the marker names in data.
#            to:  marker information (names, chr, pos) for the markers to 
#                 impute to. May be higher or lower density than from.
impute.markers = function(data, from, to) {
  # Verify that the number of row in data and from matches.
  if(nrow(data) != nrow(from)) {
    stop(paste("Number of rows in data (", nrow(data), ") must equal the",
         "number or rows in from (", nrow(from) ,")."))
  } # if(nrow(data) != nrow((from))
  # Verify that the marker names in data and from are identical.
  if(any(rownames(data) != from[,1])) {
    stop(paste("The rownames in \'data\' do not equal the marker names",
         "in column 1 of \'from\'."))
  } # if(any(rownames(data) != from[,1]))
  # If no 'to' argument was provided, jsut return the data.
  if(missing(to) || length(to) == 0) {
    return(data)
  } # if(missing(to) || length(to) == 0)
  data = as.matrix(data)
  if(!is.numeric(data)) {
    stop(paste("\'data\' must be a numeric matrix."))
  } # if(!is.numeric(data))
  # Set both marker grids to Mb.
  if(max(from[,3]) > 200) {
    from[,3] = from[,3] * 1e-6
  } # if(max(from[,3]) > 200)
  if(max(to[,3]) > 200) {
    to[,3] = to[,3] * 1e-6
  } # if(max(to[,3]) > 200)
  # Impute the values by chromosome. Use the chromosomes in the smaller set.
  chr = intersect(unique(from[,2]), unique(to[,2]))
  from = from[from[,2] %in% chr,]
  to   = to[to[,2] %in% chr,]
  newdata = matrix(0, nrow(to), ncol(data), dimnames = 
            list(to[,1], colnames(data)))
  for(c in chr) {
    from.c = from[from[,2] == c,]
    to.c   = to[to[,2] == c,]
    data.c    = data[from.c[,1],]
    newdata.c = newdata[to.c[,1],]
    
    # Interior markers.
    to.prox = outer(from.c[,3], to.c[,3], "<=")
    to.dist = outer(from.c[,3], to.c[,3], ">")
    to.prox = sapply(apply(to.prox, 2, which), max)
    to.dist = sapply(apply(to.dist, 2, which), min)
    # Trim off the ends where the array grids don't overlap.
    rng.to = which(is.finite(to.prox) & is.finite(to.dist))
    to.prox = to.prox[rng.to]
    to.dist = to.dist[rng.to]
    denom1 = from.c[to.dist,3] - from.c[to.prox,3]
    num1 = to.c[rng.to,3] - from.c[to.prox,3]
    num2 = from.c[to.dist,3] - to.c[rng.to,3]
    newdata.c[rng.to,] = (num1 * data.c[to.prox,] + num2 * data.c[to.dist,]) / 
                         denom1
    # Start of chromosome.
    # If the first marker in 'from' is > the first marker in 'to', then
    # copy the values in from[1,] to to[1,].
    to.lt.from = which(to.c[,3] < from.c[1,3])
    if(length(to.lt.from) > 0) {
      newdata.c[to.lt.from,] = matrix(newdata.c[min(rng.to),], length(to.lt.from),
                               ncol(newdata), byrow = T)
    } # if(length(to.lt.from) > 0)
    # End of chromosome.
    to.gt.from = which(to.c[,3] > from.c[nrow(from.c),3])
    if(length(to.gt.from) > 0) {
      newdata.c[to.gt.from,] = matrix(newdata.c[max(rng.to),], length(to.gt.from),
                               ncol(newdata), byrow = T)
    } # if(length(to.gt.from) > 0)
    newdata[to[,2] == c,] = newdata.c
  } # for(c)
  newdata
} # impute.markers()
