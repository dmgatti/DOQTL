################################################################################
# Given a set of founders, create all of the possible unphased genotype states 
# between them.
# Daniel Gatti
# Dan.Gatti@jax.org
# June 7, 2013
################################################################################
# Arguments: founders, character vector with founder letters.
create.genotype.states = function(founders) {
  states = outer(founders, founders, paste, sep = "")
  states = sort(states[upper.tri(states, diag = T)])
  return(list(auto = states, X = list(F = states, M = founders),
         founders = founders))
} # create.genotype.states()
