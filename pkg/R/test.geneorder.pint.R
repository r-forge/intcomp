# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.

test.geneorder.pint <- function (ge, cn.raw, cn.seg) {

  # If the ge and cn.raw are already matched
  # (are having equal resolution)
  # then always use the raw data, not segmented

  # Both raw and segmented given
  if (!is.null(cn.raw$data) && !is.null(cn.seg$data)) {
    # if raw cn resolution is same than with ge
    if (nrow(cn.raw$data) == nrow(ge$data)) {
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
