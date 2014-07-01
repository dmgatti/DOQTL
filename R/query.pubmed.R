################################################################################
# Text mine abstracts on PubMed for genes and a string (i.e. "lysosome").
# The intention is that this will be used for genes underlying QTL.
# Daniel Gatti
# Dan.Gatti@jax.org
# Oct. 22, 2013
# NOTE: Function under development.
################################################################################
# Arguments: gene.symbols: character vector containing gene symbols.
#            keywords: character vector containing words to search for in 
#                      articles involving each of the gene in gene.symbols.
query.pubmed = function(gene.symbols, keywords) {

  # Convert gene symbols to EntrezGene IDs.
  mouse.entrez = mget(gene.symbols, org.Mm.egSYMBOL2EG, ifnotfound = NA)
  mouse.entrez = lapply(mouse.entrez, function(z) { z[!is.na(z)] })
  
  # Get the human orthologs (9606 is the NCBI Homo Sapiens species code).
  message("Getting human orthologs...")
  homol = read.delim(url("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"))
  human.entrez = lapply(mouse.entrez, getHOMOLOG, targetspecies = "9606", 
                 homol = homol)
  human.entrez = lapply(human.entrez, function(z) { z[[1]][!is.na(z[[1]])] })
  human.entrez = human.entrez[sapply(human.entrez, length) > 0]
  human.entrez = lapply(human.entrez, as.character)
  
  # Get the PubMed IDs for this gene.
  mouse.pmid = lapply(mouse.entrez, mget, envir = org.Mm.egPMID, ifnotfound = NA)
  human.pmid = lapply(human.entrez, mget, envir = org.Hs.egPMID, ifnotfound = NA)
  
  # Combine all of the entrez IDs for each gene.
  pmid = as.list(rep(NA, length(mouse.entrez)))
  names(pmid) = names(mouse.entrez)
  for(i in names(pmid)) {
    pmid[[i]] = unique(c(unlist(mouse.pmid[[i]]), unlist(human.pmid[[i]])))
  } # for(i)
  
  pmid = lapply(pmid, function(z) { z[!is.na(z)] })
  message("Retrieving PubMed articles ...")
  results = as.list(gene.symbols)
  names(results) = gene.symbols
  
  for(g in 1:length(gene.symbols)) {
    message(gene.symbols[g])
    results[[g]] = vector(mode = "list", length = length(keywords))
    names(results[[g]]) = keywords
    results[[g]] = lapply(results[[g]], as.list)
    x = pubmed(pmid[[g]])
    a = xmlRoot(x)
    num.art = length(xmlChildren(a))
 
    art.list = vector("list", length = num.art)
    for(i in 1:num.art) {
      art.list[[i]] = buildPubMedAbst(a[[i]])
      found = sapply(keywords, grep, x = abstText(art.list[[i]]))
      found = which(sapply(found, length) > 0)
      if(length(found) > 0) {
        for(k in 1:length(found)) {
          results[[g]][[found[k]]] = c(results[[g]][[found[k]]], art.list[[i]])
        } # for(k)
      } # if(length(found) > 0)
    } # for(i)
  } # for(g)
  
  return(results)
  
} # query.pubmed
