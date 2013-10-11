get.machine.precision = function() {
  return(log10(.Machine$double.eps) / log10(exp(1)))
}
