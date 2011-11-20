# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.

test.geneorder.sim <- function (ge, cn, meth = "full", runname = "simtest", regs = 1, win = 1e6) {
  
  require(SIM)
  cn$data <- as.matrix(cn$data)
  ge$data <- as.matrix(ge$data)
  
  # assign absolute start positions by adding a large number in front of
  # plain location, which is then increased with chromosome index
  chr.starts <- c()
  abs.starts <- c(); n <- 1
  for (chrx in 1:22) {
    starts <- cn$info[cn$info[["chr"]] == chrx, "loc"]
    chr.starts <- c(chr.starts, starts)
    abs.starts <- c(abs.starts, n + starts)
    n <- max(abs.starts) + 1
  }
  abs.starts <- floor(abs.starts)
  chr.starts <- floor(chr.starts)

  acgh.data <- as.data.frame(list(rownames(cn$data), rownames(cn$data), cn$info[, "chr"], I(chr.starts), I(abs.starts), cn$data))
  
  colnames(acgh.data) <- c("ID", "Symbol", "CHROMOSOME", "STARTPOS",
                           "Abs.start", colnames(cn$data))
                           
  exprs.data <- as.data.frame(list(rownames(ge$data), rownames(ge$data), ge$info[, "chr"], I(chr.starts), I(abs.starts), cbind(ge$data)))
  
  colnames(exprs.data) <- c("ID", "Symbol", "CHROMOSOME", "STARTPOS", "Abs.start", colnames(ge$data))
    
  assemble.data(dep.data = acgh.data, 
                indep.data = exprs.data, 
                dep.ann = colnames(acgh.data)[1:4], 
                indep.ann = colnames(exprs.data)[1:4], 
                dep.id = "ID", 
                dep.chr = "CHROMOSOME", 
                dep.pos = "STARTPOS", 
                dep.symb = "Symbol", 
                indep.id = "ID", 
                indep.chr = "CHROMOSOME", 
                indep.pos = "STARTPOS", 
                indep.symb = "Symbol", 
                overwrite = TRUE, 
                run.name = runname)

  if (meth == "window") {
    integrated.analysis(samples = 6:(6+ncol(ge$data)-1), method = meth, run.name = runname, zscores = TRUE, input.regions = regs, window = c(win, win))
  }
  if (meth == "full") {
    integrated.analysis(samples = 6:(6+ncol(ge$data)-1), method = meth, run.name = runname, zscores = TRUE, input.regions = regs)
  }

  table.dep <- tabulate.top.dep.features(method = meth, run.name = runname, input.regions = regs, adjust.method="none")
  
  # List p-values for genes
  if (length(regs) > 1) {
    #LL's code for real data with 22 regions
    pstab <- NULL
    for (i in regs) {
      ps <- table.dep[[i]][["P.values"]]
      ids <- as.character(table.dep[[i]][["ID"]])
      pstab <- rbind(pstab, cbind(ps, ids))
    }
    
    ps <- pstab[,1]
    names(ps) <- pstab[,2]
    
    # order genes
    ids <- names(sort(ps))
  } else if (length(regs) == 1) {  
    # MS's code for simulated data
    ids <- as.character(table.dep[[regs]][["ID"]])
    # genes are already ordered 
  }
  
  # genes are already ordered 
  ids
  
}
