addLogVector = function(x) {
  # Get the largest value and save it.
  max.index = which.max(x)
  retval = x[max.index]
  # Divide all values by the largest value.
  x = x - retval
  # Set those values less than the machine precision that 
  # are too small to influence the largest value to NA.
  x[x < get.machine.precision()] = NA
  # Set the largest value to NA because we already have it in retval.
  x[max.index] = NA
  # Return the largest value * the sum of the remaining values.
  return(retval + log1p(sum(exp(x), na.rm = T)))
}

