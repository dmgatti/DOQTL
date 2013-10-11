################################################################################
# Read in the genotype.probx.txt files and write them out as *.Rdata objects.
# Arguments: prob.files: character vector with prsmth filenames.
create.Rdata.files = function(prob.files) {
  for(i in 1:length(prob.files)) {
    print(prob.files[i])
    prsmth = read.delim(prob.files[i])
    save(prsmth, file = sub("txt", "Rdata", prob.files[i]))
  } # for(i)
} # create.Rdata.files()

