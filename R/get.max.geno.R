get.max.geno = function(probs) {
  max.geno = matrix("", dim(probs)[[1]], dim(probs)[[3]], dimnames = 
             list(dimnames(probs)[[1]], dimnames(probs)[[3]]))
  founders = dimnames(probs)[[2]]
  for(s in 1:dim(probs)[[1]]) {
    if(s %% 50 == 0) print(s)
    mxp = apply(probs[s,,], 2, max)
    wmxp = apply(probs[s,,], 2, which.max)
    homo = which(mxp >= 0.75)
    max.geno[s,homo] = paste(founders[wmxp[homo]], founders[wmxp[homo]], sep = "")
    het = which(mxp < 0.75)
    rnk = apply(probs[s,,het], 2, rank)
    rnk = apply(rnk, 2, ">", 6)
    rnk = apply(rnk, 2, which)
    let = matrix("", nrow(rnk), ncol(rnk))
    let[1,] = founders[rnk[1,]]
    let[2,] = founders[rnk[2,]]
    let = apply(let, 2, sort)
    max.geno[s,het] = paste(let[1,], let[2,], sep = "")
  } # for(s)
  return(max.geno)
} # get.max.geno()
