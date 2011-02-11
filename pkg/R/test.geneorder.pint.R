test.geneorder.pint <- function (ge, cn, Labels = NULL) {

  # Note that pint does not utilize Label information, option added for compatibility
  
  message("ordering genes with pint")
  require(pint)

  # Add artificially arm name as this is required by pint
  #ge.pint <- ge; cn.pint <- cn
  #if (!"arm" %in% colnames(ge$info)) {ge.pint <- ge; ge.pint$info[["arm"]] <- rep("p", nrow(ge$info))}
  #if (!"arm" %in% colnames(cn$info)) {cn.pint <- cn; cn.pint$info[["arm"]] <- rep("p", nrow(cn$info))}
  #ge <- ge.pint
  #cn <- cn.pint
  
  # Calculate the models
  models <- screen.cgh.mrna(ge, cn)

  # Order the genes
  #ordg <- topGenes(models, n = nrow(ge$data))
  ordg <- orderGenes(models)[, "genes"]
  
  # Return ordered gene list
  ordg
}
