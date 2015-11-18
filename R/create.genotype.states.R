################################################################################
# Given a set of founders, create all of the possible unphased genotype states 
# between them.
# Daniel Gatti
# Dan.Gatti@jax.org
# June 7, 2013
################################################################################
# Arguments: founders, character vector with founder letters.
#            sampletype: one of "DO, "HS", "CC", or "HSrat"
create.genotype.states = function(founders, sampletype = c("DO", "CC", "HS", "HSrat")) {

  sampletype = match.arg(sampletype)

  states = NULL
  if(sampletype == "CC") {

    states = paste0(founders, founders)

  } else {

    states = outer(founders, founders, paste, sep = "")
    states = sort(states[upper.tri(states, diag = TRUE)])

  } # else

  return(list(auto = states, X = list(F = states, M = founders),
         founders = founders))

} # create.genotype.states()
