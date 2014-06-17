rankZ <-
function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (length(x) + 1)
  return(qnorm(x))
}
