################################################################################
# Create initial transition matrices for each SNP.  These are based on code by
# Karl Broman that uses the distance between each SNP, the DO generation and
# the pre-CC progenitors used to create the DO.
# There are 9 unique transition cases in the DO on autosomes and the femaleX.
# There are 2 on the male X. We calculate the probability of seeing a 
# recombination between each SNP for each type of transition and then place 
# the probabilities in the correct cells of the SNP-specific transition matrix.
# WARNING: THIS FUNCTION IS COMPLETELY DO SPECIFIC BECAUSE IT USES THE DO 
#          PROGENITORS AND DO GENERATION!
# Arguments: states: a list of states that are two letter combinations of the
#                    founders.
#            snp.loc: matrix with 4 columns and num.snps rows. column 1: SNPID,
#                     column 2: chr, column 3: base location, column 4: cM
#                     location.  These should only be the SNPs on the current
#                     chromosome.
#            do.gen: The DO outbreeding generation for each sample.
#            chr: one of "auto", "femaleX", "maleX", for each type of chr.
# Returns: list of transition probabilitry arrays, one for each DO generation
#          in the data set. Each list element is an array of transition
#          probabilties for each HMM state and each SNP on the current Chr.
#          Array dimensions: num.states x num.states x num.SNPs.
#          Probabilities are returned on a log scale. 
create.log.transition.matrices <-
function(states, snp.loc, do.gen, chr = c(1:19, "X"), sex = c("M", "F")) {

  chr = match.arg(chr)
  sex = match.arg(sex)

  # The original pre-CC lines used to create the DO (from K. Svenson).
  # The names are the CC generation number and the values are the number
  # of mice from that generation.
  alpha = c(21, 64, 24, 10, 5, 9, 5, 3, 3)
  alpha = alpha / sum(alpha)
  names(alpha) = c(4, 5, 6, 7, 8, 9, 10, 11, 12)

  # Get the unique DO generations.
  unique.gen = sort(unique(do.gen))
  unique.gen = unique.gen[!is.na(unique.gen)]

  # Create the return value list.
  retval = as.list(1:length(unique.gen))
  names(retval) = unique.gen

  # Create an empty list of transition probabilities.
  trans.prob = as.list(unique.gen)
  names(trans.prob) = unique.gen

  # This holds the cas ID of each type of transition probability.
  cases = matrix(0, length(states), length(states), dimnames = list(states,
                 states))

  # Females
  if(sex == "F") {
      cases = assign.autosome.femaleX.cases(states)

      for(i in 1:length(unique.gen)) {
        # Create an array of num.states x num.states x num.snps.  Each row is the
        # probability that the state in row i will transition to the state in 
        # column j.
        retval[[i]] = array(0, dim = c(length(states), length(states),
                       nrow(snp.loc) - 1), dimnames = list(states, states,
                       snp.loc[-nrow(snp.loc),1]))
      } #for(i)

  } else {
    # Males
    # Autosomes
    if(chr %in% 1:19) {
      cases = assign.autosome.femaleX.cases(states)

      for(i in 1:length(unique.gen)) {
        # Create an array of num.states x num.states x num.snps.  Each row is the
        # probability that the state in row i will transition to the state in 
        # column j.
        retval[[i]] = array(0, dim = c(length(states), length(states),
                       nrow(snp.loc) - 1), dimnames = list(states, states,
                       snp.loc[-nrow(snp.loc),1]))
      } #for(i)
    } else {

      # Male X chromosome.
      founders = sort(unique(unlist(strsplit(states, split = ""))))
      # There are only two possible changes since the male X Chr is hemizygous.
      cases = matrix(2, length(founders), length(founders))
      diag(cases) = 1

      # Create a transition probability array for each DO generation.
      for(i in 1:length(unique.gen)) {
        # Create an array of num.states x num.states x num.snps.  Each row is the
        # probability that the state in row i will transition to the state in 
        # column j.
        retval[[i]] = array(0, dim = c(length(founders), length(founders),
                      nrow(snp.loc) - 1), dimnames = list(founders, founders,
                      snp.loc[-1,1]))
      } #for(i)
    } # else
  } # else

  # Get the transition probabilities and place them in the correct locations.
  for(i in 1:length(unique.gen)) {
    trans.prob = get.trans.probs(r = diff(snp.loc[,4]) * 1e-8, do.gen = 
                 unique.gen[i], alpha, chr, sex)

    for(s in 1:(nrow(snp.loc) - 1)) {
      retval[[i]][,,s] = trans.prob[s,cases]
    } # for(s)

    retval[[i]] = log(retval[[i]])

    # This takes care of probabilities that were = 0.
    retval[[i]][is.infinite(retval[[i]])] = -.Machine$double.xmax

  } #for(i)

  return(retval)

}

