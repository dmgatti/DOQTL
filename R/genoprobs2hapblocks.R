################################################################################
# Given a DO 36-state genoprobs matrix, output a UNC-style hap file that lists
# the start and end of pseudophased haplotype blocks.
# Daniel Gatti
# dan.gatti@jax.org
# Aug. 24, 2016
################################################################################
genoprobs2hapblocks = function(probs, markers) {

  mat = get.diplotype2haplotype.matrix(get.do.states())
  probs = probs %*% mat
  markers = markers[match(rownames(probs), markers[,1]),]
  stopifnot(markers[,1] == rownames(probs))

  haps = haploprobs2hapblocks(probs, markers)

} # genoprobs2hapblocks()


haploprobs2hapblocks = function(probs, markers) {

  probs = probs[rownames(probs) %in% markers[,1],]
  markers = markers[markers[,1] %in% rownames(probs),]
  probs = probs[markers[,1],]

  # Multiply by 2 and round to get the clear haplotype calls.
  p2 = round(2 * probs)

  # Deal with the ambiguous calls (rows that don't sum to '2')
  wh = which(rowSums(p2) != 2)
  for(i in wh) {
    x = sort(probs[i,])[7:8]
    p2[i,] = 0
    p2[i,names(x)] = 1
  } # for(i)
  stopifnot(rowSums(p2) == 2)

  # Split up by chromosome.
  p2 = split(data.frame(p2), markers[,2])
  m2 = split(markers[,3], markers[,2])
  diff = p2
  haps = vector("list", length(diff))
  names(haps) = names(diff)

  # Loop through and get the haplotype boundaries on each chromosome.
  for(i in 1:length(p2)) {

    p2[[i]] = as.matrix(p2[[i]])
    rownames(p2[[i]]) = m2[[i]]
    diff[[i]] = diff(p2[[i]])
    diff[[i]] = rbind(p2[[i]][1,], diff[[i]], -p2[[i]][nrow(p2[[i]]),])
    rownames(diff[[i]])[1] = rownames(p2[[i]])[1]
    rownames(diff[[i]])[nrow(diff[[i]])] = rownames(p2[[i]])[nrow(p2[[i]])]
    diff[[i]] = diff[[i]][rowSums(diff[[i]] != 0) > 0,]

    # Check for colSums = 0, and rowSums = 0, except for first and last rows.
    stopifnot(all(colSums(diff[[i]]) == 0))
    rs = rowSums(diff[[i]])
    stopifnot(all(rs[c(-1, -length(rs))] == 0))
    stopifnot(rs[1] == 2)
    stopifnot(rs[length(rs)] == -2)
    rm(rs)

    # Create the blocks for strand 1.
    haps[[i]] = vector("list", 2)
    row = 1
    next.hap = which(diff[[i]][row ,] > 0)[1]
    start = as.numeric(rownames(diff[[i]])[row])
    diff[[i]][row, next.hap] = diff[[i]][row, next.hap] - 1

    while(!is.na(next.hap)) {

	row = min(which(diff[[i]][,next.hap] < 0 & 1:nrow(diff[[i]]) > row))
      haps[[i]][[1]] = c(haps[[i]][[1]], names(next.hap), start * 1e6 + 1, 
                       as.numeric(rownames(diff[[i]])[row]) * 1e6)
      diff[[i]][row, next.hap] = diff[[i]][row, next.hap] + 1
      next.hap = which(diff[[i]][row,] > 0)[1]
      diff[[i]][row, next.hap] = diff[[i]][row, next.hap] - 1
      start = as.numeric(rownames(diff[[i]])[row])

    } # while(!is.na(next.hap))

    if(any(is.na(haps[[i]][[1]]))) { stop(paste(i, "haps1")) }

    diff[[i]] = diff[[i]][rowSums(diff[[i]] != 0) > 0,]

    # Create the blocks for strand 2.
    next.hap = which(diff[[i]][1,] > 0)
    start = as.numeric(rownames(diff[[i]])[1])
    for(row in 2:nrow(diff[[i]])) {

      haps[[i]][[2]] = c(haps[[i]][[2]], names(next.hap), start * 1e6 + 1, 
                       as.numeric(rownames(diff[[i]])[row]) * 1e6) 
      next.hap = which(diff[[i]][row,] > 0)
      start = as.numeric(rownames(diff[[i]])[row])

    } # for(row)

    if(any(is.na(haps[[i]][[2]]))) { stop(paste(i, "haps2")) }

  } # for(i)

  # Resort the chromosomes so they're in numeric order.
  haps = haps[order(as.numeric(names(haps)))]

  return(haps)

} # haploprobs2hapblocks()


write.unc.hap.file = function(haps, filename, title) {

  outf = file(description = filename, open = "w")
  writeLines(text = c("set,ppc,20"), con = outf)
  writeLines(text = paste0("title,top,", title), con = outf)
  writeLines(text = "title,bottom,\"Version 2.2\"", con = outf)
  writeLines(text = "grid,2.0", con = outf)
  writeLines(text = "set,overlap,sticky", con = outf)

  for(i in 1:length(haps)) {

    haps[[i]][[1]] = paste(haps[[i]][[1]], collapse = ",")
    writeLines(text = paste0("chr", names(haps)[i], ",,", haps[[i]][[1]]),
               con = outf)
    writeLines(text = paste0("gap0,\"", names(haps)[i], "\""), con = outf)
    haps[[i]][[2]] = paste(haps[[i]][[2]], collapse = ",")
    writeLines(text = paste0("chr", names(haps)[i], ",,", haps[[i]][[2]]),
               con = outf)
    writeLines(text = "gap", con = outf)

  } # for(i)

  close(outf)

} # write.unc.hap.file()

