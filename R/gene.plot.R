################################################################################
# Given a genomic range, plot the genes in the interval.
# Daniel Gatti
# Dan.Gatti@jax.org
# Feb. 13, 2013
################################################################################
# Arguments: mgi.file: data.frame, containing genes (or other features) as 
#                      returned by get.mgi.features().
#            rect.col: color, the color of the gene rectangles.
#            text.col: color, the color of the gene names.
#            ...: other arguments to be passed to plot.
# Returns: data.frame with gene symbol, gene start, gene end, text start, 
#          text end and row.
gene.plot = function(mgi, rect.col = "grey30", text.col = "black", ...) {
  # If we have no genes, just plot an empty frame and return.
  if(is.null(mgi) || length(mgi) == 0) {
    plot(0, 0, col = 0, xlab = "", xaxs = "i", ylab = "", yaxt = "n", ...)
    return()
  } # if(is.null(mgi) || nrow(mgi) == 0)
  
  previous.cex = par("cex")
  call = match.call(expand.dots = TRUE)
  mgi$Name = as.character(mgi$Name)
  
  # Convert the mgi gene start and stop locations to Mb.
  if(all(mgi$start > 200)) {
     mgi$start = mgi$start * 1e-6
  } # if(all(mgi$start > 200))
  if(all(mgi$stop > 200)) {
     mgi$stop = mgi$stop * 1e-6
  } # if(all(mgi$stop > 200))
  # Place names on the colors so that we can keep track of which
  # colors are associated with which genes. Also, expand them if they are
  # shorter than nrow(mgi).
  if(length(rect.col) < nrow(mgi)) {
    rect.col = rep(rect.col, ceiling(nrow(mgi) / length(rect.col)))[1:nrow(mgi)]
  } # if(length(rect.col) < nrow(mgi))
  names(rect.col) = mgi$Name
  if(length(text.col) < nrow(mgi)) {
    text.col = rep(text.col, ceiling(nrow(mgi) / length(text.col)))[1:nrow(mgi)]
  } # if(length(text.col) < nrow(mgi))
  names(text.col) = mgi$Name
  
  # If we have xlim in the arguments, then make the plot using the user
  # defined xlim and subset the mgi data to include only genes within the
  # plot limits.
  if(any(names(call) == "xlim")) {
    xlim = eval(call$xlim, envir = parent.frame())
    plot(0, 0, col = 0, xlab = "", xaxs = "i", ylab = "", yaxt = "n", ...)
    mgi = mgi[mgi$stop >= xlim[1] & mgi$start <= xlim[2],]
    # If we have no genes to plot, just return.
    if(nrow(mgi) == 0) {
      return()
    } # if(nrow(mgi) == 0)
  } else {
    plot(0, 0, col = 0, xlim = c(min(mgi$start), max(mgi$stop)),
         xlab = "", ylab = "", yaxt = "n", ...)
  } # else
  mtext(side = 1, line = 2, text = paste("Chr", mgi$seqid[1], "(Mb)"))
  # Line the genes up sequentially in columns.
  # Locs holds the gene symbol, the gene start and end, the text start and end,
  # as well as the row to plot on. Use strwidth("W") as the space after a gene.
  locs = data.frame(name = mgi$Name, gstart = mgi$start,
         gend = mgi$stop, tstart = mgi$stop + strwidth("i"), 
         tend = mgi$stop + strwidth("W") + sapply(mgi$Name, strwidth),
         stringsAsFactors = F)
  locs = locs[order(locs$gstart),]
  par(lend = 2)
  usr = par("usr")
  nrows = 1 # Number of rows in the current plot.
  row = 10  # Number of rows that we need.
  ymin = 1  # Lowest y-value for the gene closest to the bottom of the plot.
  iter = 0  # Number of iterations.
  # We need at least enough rows in the plot to fit the data.
  while((nrows < row | (usr[4] - ymin) / diff(usr[3:4]) < 0.7) & iter < 20) { 
    last.strht = strheight("I") 
    retval = get.gene.locations(locs, usr)
    ymin  = retval$newloc$bottom[nrow(retval$newloc)]
    nrows = retval$nrows
    row   = retval$row
    # If we need more rows in the plot, then decrement cex until strheight("I") decreases.
    if(row > nrows) {
      while(strheight("I") == last.strht) {
        par(cex = par("cex") * 0.98)
      } # while(strheight("I") == last.strht)
    } else if (row < nrows) {
      while(strheight("I") == last.strht) {
        par(cex = par("cex") * 1.02)
      } # while(strheight("I") == last.strht)
    } # else if (row < nrows)
    iter = iter + 1
  } # while((nrows < row | (usr[4] - ymin) ...
  # If we have any horizontal collisions, shrink the font size slightly.
  if(any((retval$newloc$textx + strwidth(retval$newloc$name))[-1] -
         retval$newloc$left[-nrow(retval$newloc)] < 0)) {
    last.strht = strheight("I") 
    while(strheight("I") == last.strht) {
      par(cex = par("cex") * 0.99)
    } # while(strheight("I") == last.strht)
    retval = get.gene.locations(locs, usr)
  } # if(any((retval$newloc$textx  ...
  # Plot the genes.
  rect(retval$newloc$left, retval$newloc$bottom, retval$newloc$right, 
       retval$newloc$top, col = rect.col[retval$newloc$name],
       border = rect.col[retval$newloc$name])  
  text(retval$newloc$textx, retval$newloc$texty, retval$newloc$name, 
       col = text.col[retval$newloc$name], adj = c(0,0))
  
  par(cex = previous.cex)
  return(locs)
} # gene.plot()
# Helper function to set gene locations on plot.
get.gene.locations = function(locs, usr) {
    boxheight = 1.1 * strheight("I")
    offset = 0.15 * boxheight
    rowheight = boxheight + offset
    nrows = floor(diff(usr[3:4]) / rowheight)
    row = 1
    x = min(locs$gstart)
    tmp = locs # tmp is a sacrificial data frame from which we will remove plotted genes.
    # Try to fill in the genes without collisions.
    # newloc holds the positions of the gene rectangles and names in user coordiantes.
    newloc = data.frame(name = rep("", nrow(locs)), left = rep(0, nrow(locs)), 
             bottom = rep(0, nrow(locs)), right = rep(0, nrow(locs)), 
             top = rep(0, nrow(locs)), textx = rep(0, nrow(locs)), 
             texty = rep(0, nrow(locs)), stringsAsFactors = F)
    i = 1    # Gene counter.
    while(i <= nrow(locs)) {
      # Find the gene with the minimum start position past the current X position.
      wh = which(tmp$gstart >= x)
      # If we found one, then add it to our list of positions (newloc).
      if(length(wh) > 0) {
        # idx indexes the current gene row in tmp.
        idx = min(wh)
        newloc$name[i] = tmp$name[idx]
        newloc$left[i] = tmp$gstart[idx]
        newloc$bottom[i] = usr[4] - row * rowheight + offset
        newloc$right[i] = tmp$gend[idx]
        newloc$top[i] = usr[4] - (row - 1) * rowheight - offset
        newloc$textx[i] = tmp$tstart[idx]
        newloc$texty[i] = usr[4] - row * rowheight + offset
        x = newloc$textx[i] + strwidth(newloc$name[i]) + strwidth("i")
        # Remove this gene from tmp.
        tmp = tmp[-idx,]
        # If our X position has moved past the right end of the plot, advance to the next
        # row and reset the X position to the left edge of the plot (in user coordinates).
        if(x > usr[2]) {
          row = row + 1
          x = min(locs$gstart)
        } # if(x > usr[2])
        # Increment the gene counter.
        i = i + 1
      } else {
        # If we didn't find a gene past out current X position, advance to the next
        # row and reset the X position to the left edge of the plot (in user coordinates).
        row = row + 1
        x = min(locs$gstart)
      } # else
    } # while(i <= nrow(locs))
    return(list(nrows = nrows, row = row, newloc = newloc))
} # get.gene.locations()
