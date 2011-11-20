# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.

list.methods <- function (Labels = c(1,2)) {
  res <- test.geneorder.pipeline(ge = list(data = matrix(rnorm(9),3,3)),
                                 cn.raw = list(data = matrix(rnorm(9),3,3)), 
          cancerGenes = NULL, evaluate = FALSE, cn.default = "raw")
  res$available.methods
}