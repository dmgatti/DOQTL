get.max.geno = function(probs) {

  max.geno = matrix("", nrow(probs), dim(probs)[[3]], dimnames = 
             list(rownames(probs), dimnames(probs)[[3]]))
  founders = dimnames(probs)[[2]]

  if(dim(probs)[3] == 1) {
  
    mxp  = apply(probs[,,1], 1, max)
    wmxp = apply(probs[,,1], 1, which.max)
    homo = which(mxp >= 0.75)
    max.geno[homo,1] = paste(founders[wmxp[homo]], founders[wmxp[homo]], sep = "")
    het = which(mxp < 0.75)
    rnk = apply(probs[het,,], 1, rank)
    rnk = apply(rnk, 2, ">", 6)
    rnk = apply(rnk, 2, which)
    let = matrix("", nrow(rnk), ncol(rnk))
    let[1,] = founders[rnk[1,]]
    let[2,] = founders[rnk[2,]]
    let = apply(let, 2, sort)
    max.geno[het,1] = paste(let[1,], let[2,], sep = "")

  } else {
  
    # Loop through each sample.
    for(s in 1:nrow(probs)) {
  
      if(s %% 50 == 0) print(s)
      mxp  = apply(probs[s,,,drop = FALSE], 3, max)
      wmxp = apply(probs[s,,,drop = FALSE], 3, which.max)
      homo = which(mxp >= 0.75)
      if(length(homo) > 0) {
        max.geno[s, homo] = paste(founders[wmxp[homo]], founders[wmxp[homo]], sep = "")
      } # if(length(homo) > 0)
      het = which(mxp < 0.75)
      if(length(het) > 0) {
        rnk = apply(probs[s,,het,drop = FALSE], 3, rank)
        rnk = apply(rnk, 2, ">", 6)
        rnk = apply(rnk, 2, which)
        let = matrix("", nrow(rnk), ncol(rnk))
        let[1,] = founders[rnk[1,]]
        let[2,] = founders[rnk[2,]]
        let = apply(let, 2, sort)
        max.geno[s, het] = paste(let[1,], let[2,], sep = "")
      } # if(length(het) > 0)
    
    } # for(s)

  } # else
  
  return(max.geno)
  
} # get.max.geno()
