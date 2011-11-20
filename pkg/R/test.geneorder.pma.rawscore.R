# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.

test.geneorder.pma.rawscore <- function (ge, cn) {

  # Note that this method does not utilize Label information, option added for compatibility

  # Use original scores by PMA to order the genes

  ## use type="ordered"
  ## since both gene expression and copy number data sets are
  ## ordered by their chromosomal location

  require(PMA)
  
  tab <- get.pma.table(cn, ge, typex = "ordered", typez = "ordered")
        
  #############################################
	       
  # Order genes by PMA score. Higher scores first etc.
  #tab.ordered <- esort(tab, -score)
  scores <- tab[, "score"]
  tab.ordered <- tab[order(scores, decreasing = TRUE),]
		     
  # order genes
  ordg <- tab.ordered[, "probeid"]

  ordg
}


