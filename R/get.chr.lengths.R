#############################################################################
# Get the mouse chromosome lengths.
get.chr.lengths = function() {
  # Get the chromosome lengths.
  chrlen = org.Mm.egCHRLENGTHS
  chrlen = chrlen[names(chrlen) %in% c(1:19, "X", "Y", "M")]
  names(chrlen) = sub("X", "20", names(chrlen))
  names(chrlen) = sub("Y", "21", names(chrlen))
  names(chrlen) = sub("M", "22", names(chrlen))
  chrlen = chrlen[order(as.numeric(names(chrlen)))]
  names(chrlen) = sub("20", "X", names(chrlen))
  names(chrlen) = sub("21", "Y", names(chrlen))
  names(chrlen) = sub("22", "M", names(chrlen))
  chrlen = chrlen * 1e-6
  return(chrlen)
} # get.chr.lengths()
