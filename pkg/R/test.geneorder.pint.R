test.geneorder.pint <- function (ge, cn.raw, cn.seg, Labels = NULL) {

  # If the ge and raw cn have same resolution
  # then use the raw data, not segmented: pint will 
  # automatically detect the dependent variation.

  # Both raw and segmented given
  if (!is.null(cn.raw$data) && !is.null(cn.seg$data)) {
    # if raw cn resolution is same than with ge
    if (nrow(cn.raw$data)/nrow(ge$data) == 1) {
      cn <- cn.raw    
    } else {
      cn <- cn.seg
    }
  } else if (!is.null(cn.seg)) {
    cn <- cn.seg
  } else if (!is.null(cn.raw)) {
    cn <- cn.raw
  }
  
  # Note that pint does not utilize Label information, option added for compatibility
  require(pint)
  
  # Calculate the models
  models <- screen.cgh.mrna(ge, cn)

  # Order the genes
  ordg <- orderGenes(models)[, "genes"]
  
  # Return ordered gene list
  ordg
  
}
