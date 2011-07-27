test.geneorder.pint <- function (ge, cn, Labels = NULL) {

  # Note that pint does not utilize Label information, option added for compatibility
  
  message("ordering genes with pint")
  require(pint)
  
  # Calculate the models
  models <- screen.cgh.mrna(ge, cn, match.probes = FALSE) # matched already

  # Order the genes
  #print("Order GENES")
  #save(models, file = "~/tmp/tmp.RData")
  ordg <- orderGenes(models)[, "genes"]
  
  # Return ordered gene list
  ordg
}
