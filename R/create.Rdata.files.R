create.Rdata.files = function(prob.files) {
  for(i in 1:length(prob.files)) {
    print(prob.files[i])
    prsmth = read.delim(prob.files[i])
    prsmth = exp(as.matrix(prsmth))
    class(prsmth) = c("genoprobs", class(prsmth))
    save(prsmth, file = sub("txt", "Rdata", prob.files[i]))
  } # for(i)
} # create.Rdata.files()
