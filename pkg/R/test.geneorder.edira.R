# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.

test.geneorder.edira <- function(ge, cn, references){
  require(edira)
  # calculate edira
  edi <- edira_ratios(ge, cn, references)
  # order genes
  names(edi$test) <- rownames(ge$data)
  ordg <- names(sort(edi$test))  
  ordg
}
