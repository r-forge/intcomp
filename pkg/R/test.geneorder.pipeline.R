# ge            list object contain a data matrix with gene expression data and a matrix with annotations for gene expression data
# cn            object contain a data matrix with copy number data and a matrix with annotations for copy number data
# Labels        labels to distinguish between patient (1) and reference (0) samples
# cancerGenes   character vector containing the labels of the cancer genes
# nperm         number of permutations to perform
# input         are the data to be analysed real or simulated data?
# version       (applies to some methods only) "normal" performs gene-sampling, while "approximate" performs sample-sampling
# methods       which methods should be performed? at least one of c("edira","DRI.cp","DRI.cs","DRI.ct","SIM.full","SIM.window","intcngean","PMA","PMA.raw","pint","PREDA","DRI.ss","DRI.srank","DRI.sraw")

test.geneorder.pipeline <- function (ge, cn.raw, cn.seg=NULL,
cn.call=NULL, cghCall = NULL, ge.norm = NULL, cn.norm = NULL, Labels=NULL,
cancerGenes, nperm = 1e2, input="real", version = "normal", methods =
NULL, chromosomes = as.character(1:22)) {

  # Determine default cn to be used in all methods unless specified otherwise
  # (see intCNGEan and CNAmet)
  # By defalt use cn.raw; if not given use cn.seg; if not given, use cn.call
  if (!is.null(cn.raw)) {
    cn <- cn.raw
  } else if (is.null(cn.raw) && !is.null(cn.seg)) {
    cn <- cn.seg
  } else if (is.null(cn.raw) && is.null(cn.seg) && !is.null(cn.call)) {
    cn <- cn.call
  }
  
  ###########################
      
  ge2 <- list(data=cbind(ge$data,ge.norm$data), info=ge$info) 
  cn2 <- list(data=cbind(cn$data,cn.norm$data), info=cn$info)
    
  # If no labels given, use the same Label for all samples
  if (is.null(Labels)) { Labels <- rep(1, ncol(ge$data)) }
  
  auc <- list()

  if (!is.null(methods) && ("edira" %in% methods)) {
      message("edira")
      ordg <- test.geneorder.edira(ge, cn, Labels)
      auc[["edira"]] <- roc.auc(ordg, cancerGenes)
  }

  if (!is.null(methods) && ("DRI.cp" %in% methods)) {
    message("DRI.cp")
    ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="pearson", version=version)
    auc[["DRI.cp"]] <- roc.auc(ordg, cancerGenes)
  }
  
  if (!is.null(methods) && ("DRI.cs" %in% methods)) {
    message("DRI.cs")
    ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="spearman", version=version)
    auc[["DRI.cs"]] <- roc.auc(ordg, cancerGenes)
  }

  if (!is.null(methods) && ("DRI.ct" %in% methods)) {
      message("DRI.ct")
      ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="ttest", version=version)
      auc[["DRI.ct"]] <- roc.auc(ordg, cancerGenes)
  }
 
  if (!is.null(methods) && ("SIM.full" %in% methods)) {
    message("SIM.full")

    if(input == "real"){
      ordg <- test.geneorder.sim(ge, cn, meth = "full", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = 1:22)
    }
    if(input == "simulations.equal.dimensions" || input == "simulations.unequal.dimensions"){
      ordg <- test.geneorder.sim(ge, cn, meth = "full", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = 1)
    }
    auc[["SIM.full"]] <- roc.auc(ordg, cancerGenes)
  }
  
  if (!is.null(methods) && ("SIM.window" %in% methods)) {
    message("SIM.window")
    
    if(input == "real"){
      ordg <- test.geneorder.sim(ge, cn, meth = "window", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = 1:22, win = 1e6)
    }
    if(input == "simulations.equal.dimensions" || input == "simulations.unequal.dimensions"){
      ordg <- test.geneorder.sim(ge, cn, meth = "window", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = 1, win = 1e6)
    }
    auc[["SIM.window"]] <- roc.auc(ordg, cancerGenes)
  }

  if (!is.null(methods) && ("CNAmet" %in% methods)) {
     message("CNAmet")
     ordg <- test.geneorder.CNAmet(ge, cn=cn.call, Labels=NULL, nperm)
     auc[["CNAmet"]] <- roc.auc(ordg, cancerGenes)
  }
  
  if (!is.null(methods) && ("intcngean" %in% methods)) {
    message("intCNGEan")
    if(input == "real"){
      ordg <- test.geneorder.intcngean(ge, cn=cghCall, meth = "wmw", analysis.type = "univariate", 
                                       nperm = nperm, pth = 0.1)
    }                                 

    if(input == "simulations.equal.dimensions"){
      ordg <- test.geneorder.intcngean(ge, cn=cn.call, 
                                       meth="wmw",
                                       analysis.type="univariate",
                                       nperm = nperm, pth = 0.1, callprobs=sim$callprobs)
    }
    if(input == "simulations.unequal.dimensions"){
      ordg <- test.geneorder.intcngean(ge, cn=cn.call, 
                                       meth="wmw",
                                       analysis.type="univariate", 
                                       nperm = nperm, pth = 0.1, callprobs=sim$callprobs)
    }
    auc[["intCNGEan.wmw.univariate"]] <- roc.auc(ordg, cancerGenes)
  }

  if (!is.null(methods) && ("PMA" %in% methods)) {
    message("PMA")
    ordg <- test.geneorder.pma(ge, cn, Labels, nperm)
    auc[["PMA"]]  <- roc.auc(ordg, cancerGenes)
  }
  
  if (!is.null(methods) && ("PMA.raw" %in% methods)) {
    message("PMA.raw")
    ordg <- test.geneorder.pma.rawscore(ge, cn, Labels)
    auc[["PMA.raw"]]  <- roc.auc(ordg, cancerGenes)
  }
  
  if (!is.null(methods) && ("pint" %in% methods)) {
    message("pint")
    ordg <- test.geneorder.pint(ge, cn, Labels)
    auc[["pint"]] <- roc.auc(ordg, cancerGenes)
  }

  if (!is.null(methods) && ("jrivas" %in% methods)) {
    message("jrivas")
    ordg <- test.geneorder.jrivas(ge, cn, Labels)
    auc[["jrivas"]] <- roc.auc(ordg, cancerGenes)              
  }                                      
	      
  if (!is.null(methods) && ("PREDA" %in% methods)) {
    message("PREDA")
    ordg <- test.geneorder.preda(ge, cn, Labels, nperm=nperm, cancerGenes=cancerGenes,
        ge.qval.threshold=0.05, cn.qval.threshold=0.01, smoothMethod="lokern_scaledBandwidth_repeated",
        ge.smoothStatistic.threshold.up=0.5, ge.smoothStatistic.threshold.down=-0.5,
        cn.smoothStatistic.threshold.gain=0.1, cn.smoothStatistic.threshold.loss=-0.1, correction.method="fdr",
        chromosomes=unique(ge$info$chr))
    best_case <- roc.auc(ordg$best_case_order, cancerGenes)
    worst_case <- roc.auc(ordg$worst_case_order, cancerGenes)
    auc[["preda"]] <- mean(c(best_case,worst_case))
  }
  
  # Run methods below only when we have two-group comparison setup
  if (length(unique(Labels)) == 2) {

    # DRI-SAM
    if (!is.null(methods) && ("DRI.ss" %in% methods)) {
      message("DRI.ss")
      ordg <- test.geneorder.dri.sam(ge=ge2, cn=cn2, Labels, nperm=nperm, transform.type="standardize", version=version)
      auc[["DRI.ss"]] <- roc.auc(ordg, cancerGenes)
    }

    if (!is.null(methods) && ("DRI.srank" %in% methods)) {
      message("DRI.srank")
      ordg <- test.geneorder.dri.sam(ge=ge2, cn=cn2, Labels, nperm=nperm, transform.type="rank", version=version)
      auc[["DRI.srank"]] <- roc.auc(ordg, cancerGenes)
    }
    
    if (!is.null(methods) && ("DRI.sraw" %in% methods)) {
      message("DRI.sraw")
      ordg <- test.geneorder.dri.sam(ge=ge2, cn=cn2, Labels, nperm=nperm, transform.type="raw", version=version)
      auc[["DRI.sraw"]] <- roc.auc(ordg, cancerGenes)
    }

  }

  return(auc)
  
}
