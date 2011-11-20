# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.


test.geneorder.CNAmet <- function (ge, cn, nperm) {

  # NOTE: assumes called data (cn.call) in input!

  require(CNAmet)
  
  # Get amplification p-values
  cghmat <- (cn$data > 0) - 0
  results <- CNAmet(ge$data, cghmat, methylMatrix = NULL, perms = nperm, gainData = TRUE)
  pv.amp <- results[, "CWPvalue"]
   
  # Get deletion p-values
  cghmat <- (cn$data < 0) - 0
  results <- CNAmet(ge$data, cghmat, methylMatrix = NULL, perms = nperm, gainData = FALSE)
  pv.del <- results[, "CWPvalue"]
  
  # Combine amplification and deletion scores
  # If one of them is NA then use the one that has p-value
  # If both have pvals then use the smaller p-value
  # If both are NA then use NA
  pv <- c()
  for (i in names(pv.del)) {
    print(i)
    a <- pv.amp[[i]]
    d <- pv.del[[i]]
    if (is.na(a) && is.na(d)) {pv[[i]] <- NA}
    if (is.na(a) && !is.na(d)) {pv[[i]] <- d}
    if (!is.na(a) && is.na(d)) {pv[[i]] <- a}
    if (!is.na(a) && !is.na(d)) {pv[[i]] <- min(a,d)}    
  }

  # Sort by p-values
  # If both are NA add these to the end of the list with random order (about 30% in Hyman data for instance)
  print(" Order the genes")
  ordg <- c(names(sort(pv)), sample(names(which(is.na(pv)))))

  ordg
}

  
