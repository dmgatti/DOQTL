################################################################################
# Given the states, assign the 9 unique transition cases to each cell and return
# a length(states) x length(states) matrix with the case locations.
# Arguments:
# states: character vector with the 36 possible DO states.
# Returns:
# cases: numeric matrix with dimensions length(states) x length(states) that
#        contains values between 1:9 indicating which type of transition occurs
#        in each cell.
assign.autosome.femaleX.cases <- function(states) {

  # Assign transition case IDs to each cell in the transition probability matrix.
  # Copied from Karl's DOstep.c functions.
  spl = matrix(unlist(strsplit(states, split = "")), nrow = 2)
  founders = sort(unique(as.vector(spl)))
  cases = matrix(0, length(states), length(states))
  for(i in 1:length(states)) {
    for(j in 1:length(states)) {
      if(spl[1,i] == spl[2,i]) { 
        if(spl[1,j] == spl[2,j]) {
          if(spl[1,i] == spl[1,j]) {
            # AA -> AA
            cases[i,j] = 1
          } else {
            # AA -> BB
            cases[i,j] = 2
          } # else
        } else {
          if(spl[1,i] == spl[1,j] || spl[1,i] == spl[2,j]) {
            # AA -> AB
            cases[i,j] = 3
          } else {
            # AA -> BC
            cases[i,j] = 4
          } # else
        } # else
      } # else
      else {
        # AB
        if(spl[1,j] == spl[2,j]) {
          if(spl[1,i] == spl[1,j] || spl[2,i] == spl[1,j]) {
            # AB -> AA
            cases[i,j] = 5
          } else {
            # AB -> CC
            cases[i,j] = 6
          } # else
        } # else
        else {
          if((spl[1,i] == spl[1,j] && spl[2,i] == spl[2,j]) ||
             (spl[1,i] == spl[2,j] && spl[2,i] == spl[1,j])) {
            # AB -> AB
            cases[i,j] = 7
          } else if(spl[1,i] == spl[1,j] || spl[1,i] == spl[2,j] ||
                  spl[2,i] == spl[1,j] || spl[2,i] == spl[2,j]) {
            # AB -> AC
            cases[i,j] = 8
          } else { 
            # AB -> CD
            cases[i,j] = 9
          } # else
        } # else
      } # else
    } # for(j)
  } # for(i)

  return(cases)
}

