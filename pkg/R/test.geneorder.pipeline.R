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
NULL, chromosomes = as.character(1:22), callprobs, references, evaluate = TRUE) {

  # Determine default cn to be used in all methods unless specified otherwise
  # (see intCNGEan and CNAmet)
  # By defalt use cn.seg; if not given use cn.raw; if not given, use cn.call
  if (!is.null(cn.seg)) {
    cn <- cn.seg
  } else if (is.null(cn.seg) && !is.null(cn.raw)) {
    cn <- cn.raw
  } else if (is.null(cn.seg) && is.null(cn.raw) && !is.null(cn.call)) {
    cn <- cn.call
  }
  
  ###########################
  if(!is.null(ge.norm$data)) { ge2 <- list(data=cbind(ge$data,ge.norm$data), info=ge$info) }
  if(!is.null(cn.norm$data)) { cn2 <- list(data=cbind(cn$data,cn.norm$data), info=cn$info) }
    
  # If no labels given, use the same Label for all samples
  if (is.null(Labels)) { Labels <- rep(1, ncol(ge$data)) }
  
  roc <- list()
  runtime <- list()
  ordered.genes <- list()
  available.methods <- c()

  available.methods <- c(available.methods, "edira")
  if (!is.null(methods) && ("edira" %in% methods) && evaluate) {
      message("edira")
      start.time <- Sys.time()      
      ordg <- test.geneorder.edira(ge, cn, Labels, references)
      end.time <- Sys.time()
      runtime[["edira"]] <- as.numeric(difftime(end.time, start.time, units='mins'))

      roc[["edira"]] <- roc.auc2(ordg, cancerGenes)
      ordered.genes[["edira"]] <- ordg
  }

  available.methods <- c(available.methods, "DRI.cp")
  if (!is.null(methods) && ("DRI.cp" %in% methods) && evaluate) {
    message("DRI.cp")
    start.time <- Sys.time()
    ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="pearson", version=version)
    end.time <- Sys.time()
    runtime[["DRI.cp"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["DRI.cp"]] <- roc.auc2(ordg, cancerGenes)
    ordered.genes[["DRI.cp"]] <- ordg
  }
  
  available.methods <- c(available.methods, "DRI.cs")  
  if (!is.null(methods) && ("DRI.cs" %in% methods) && evaluate) {
    message("DRI.cs")
    start.time <- Sys.time()    
    ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="spearman", version=version)
    end.time <- Sys.time()
    runtime[["DRI.cs"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["DRI.cs"]] <- roc.auc2(ordg, cancerGenes)
    ordered.genes[["DRI.cs"]] <- ordg
  }
  
  available.methods <- c(available.methods, "DRI.ct")
  if (!is.null(methods) && ("DRI.ct" %in% methods) && evaluate) {
      message("DRI.ct")
      start.time <- Sys.time()      
      ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="ttest", version=version)
      end.time <- Sys.time()
      runtime[["DRI.ct"]] <- as.numeric(difftime(end.time, start.time, units='mins'))      
      roc[["DRI.ct"]] <- roc.auc2(ordg, cancerGenes)
      ordered.genes[["DRI.ct"]] <- ordg      
  }
 
  available.methods <- c(available.methods, "SIM.full")
  if (!is.null(methods) && ("SIM.full" %in% methods) && evaluate) {
    message("SIM.full")
    start.time <- Sys.time()
    if(input == "real"){
      ordg <- test.geneorder.sim(ge, cn, meth = "full", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = 1:22)
    }
    if(input == "simulations.equal.dimensions" || input == "simulations.unequal.dimensions"){
      ordg <- test.geneorder.sim(ge, cn, meth = "full", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = 1)
    }
    end.time <- Sys.time()
    runtime[["SIM.full"]] <- as.numeric(difftime(end.time, start.time, units='mins'))
    roc[["SIM.full"]] <- roc.auc2(ordg, cancerGenes)
    ordered.genes[["SIM.full"]] <- ordg    
  }
  
  available.methods <- c(available.methods, "SIM.window")  
  if (!is.null(methods) && ("SIM.window" %in% methods) && evaluate) {
    message("SIM.window")
    start.time <- Sys.time()    
    if(input == "real"){
      ordg <- test.geneorder.sim(ge, cn, meth = "window", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = 1:22, win = 1e6)
    }
    if(input == "simulations.equal.dimensions" || input == "simulations.unequal.dimensions"){
      ordg <- test.geneorder.sim(ge, cn, meth = "window", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = 1, win = 1e6)
    }
    end.time <- Sys.time()
    runtime[["SIM.window"]] <- as.numeric(difftime(end.time, start.time, units='mins'))
    roc[["SIM.window"]] <- roc.auc2(ordg, cancerGenes)
    ordered.genes[["SIM.window"]] <- ordg
  }

  available.methods <- c(available.methods, "CNAmet")
  if (!is.null(methods) && ("CNAmet" %in% methods) && evaluate) {
     message("CNAmet")
     start.time <- Sys.time()     
     ordg <- test.geneorder.CNAmet(ge, cn=cn.call, Labels=NULL, nperm)
     end.time <- Sys.time()
     runtime[["CNAmet"]] <- as.numeric(difftime(end.time, start.time, units='mins'))     
     roc[["CNAmet"]] <- roc.auc2(ordg, cancerGenes)
     ordered.genes[["CNAmet"]] <- ordg
  }
  
  available.methods <- c(available.methods, "intcngean")  
  if (!is.null(methods) && ("intcngean" %in% methods) && evaluate) {
    message("intCNGEan")
    start.time <- Sys.time()    
    if(input == "real"){
      ordg <- test.geneorder.intcngean(ge=ge,
                                       cghCall=cghCall,
                                       meth = "wmw",
                                       analysis.type = "univariate",
                                       nperm = nperm,
                                       pth = 0.1,
                                       match = TRUE)
    }                                 

    if(input == "simulations.equal.dimensions"){
      ordg <- test.geneorder.intcngean(ge=ge, cghCall=cghCall, 
                                       meth="wmw",
                                       analysis.type="univariate",
                                       nperm = nperm, pth = 0.1, callprobs=callprobs, match=FALSE)
    }
    if(input == "simulations.unequal.dimensions"){
      ordg <- test.geneorder.intcngean(ge=ge, cghCall=cghCall, 
                                       meth="wmw",
                                       analysis.type="univariate", 
                                       nperm = nperm, pth = 0.1, callprobs=NULL, match=TRUE)
    
    
    }
    end.time <- Sys.time()
    runtime[["intCNGEan.wmw.univariate"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["intCNGEan.wmw.univariate"]] <- roc.auc2(ordg, cancerGenes)
    ordered.genes[["intCNGEan.wmw.univariate"]] <- ordg
  }

  #if (!is.null(methods) && ("PMA" %in% methods)) {
  #  # NOTE: PMA.raw is the original PMA method. This PMA function contains
  #  # additional permutation test to calculate significances
  #  message("PMA")
  #  start.time <- Sys.time()    
  #  ordg <- test.geneorder.pma(ge, cn, Labels, nperm)
  #  end.time <- Sys.time()
  #  runtime[["PMA"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
  #  roc[["PMA"]]  <- roc.auc2(ordg, cancerGenes)
  #  ordered.genes[["PMA"]] <- ordg    
  #}
  
  available.methods <- c(available.methods, "PMA.raw")  
  if (!is.null(methods) && ("PMA.raw" %in% methods) && evaluate) {
    message("PMA.raw")
    start.time <- Sys.time()    
    ordg <- test.geneorder.pma.rawscore(ge, cn, Labels)
    end.time <- Sys.time()
    runtime[["PMA.raw"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["PMA.raw"]]  <- roc.auc2(ordg, cancerGenes)
    ordered.genes[["PMA.raw"]] <- ordg
  }
  
  available.methods <- c(available.methods, "pint")  
  if (!is.null(methods) && ("pint" %in% methods) && evaluate) {

    message("pint")
  
    start.time <- Sys.time()    
    ordg <- test.geneorder.pint(ge, cn.raw, cn.seg, Labels)
    end.time <- Sys.time()
    runtime[["pint"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["pint"]] <- roc.auc2(ordg, cancerGenes)
    ordered.genes[["pint"]] <- ordg    
  }

  available.methods <- c(available.methods, "OrtizEstevez")
  if (!is.null(methods) && ("OrtizEstevez" %in% methods) && evaluate) {
    message("OrtizEstevez")
    start.time <- Sys.time()    
    ordg <- test.geneorder.OrtizEstevez(ge, cn, Labels)
    end.time <- Sys.time()
    runtime[["OrtizEstevez"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["OrtizEstevez"]] <- roc.auc2(ordg, cancerGenes)
    ordered.genes[["OrtizEstevez"]] <- ordg        
  }                                      
          
  available.methods <- c(available.methods, "PREDA")	  
  if (!is.null(methods) && ("PREDA" %in% methods) && evaluate) {
    message("PREDA")
    start.time <- Sys.time()    
    ordg <- test.geneorder.preda(ge, cn, Labels, nperm=nperm, cancerGenes=cancerGenes,
        ge.qval.threshold=0.05, cn.qval.threshold=0.01, smoothMethod="spline",
        ge.smoothStatistic.threshold.up=0.5, ge.smoothStatistic.threshold.down=-0.5,
        cn.smoothStatistic.threshold.gain=0.1, cn.smoothStatistic.threshold.loss=-0.1, correction.method="fdr",
        chromosomes=unique(ge$info$chr))
    end.time <- Sys.time()    
    runtime[["preda"]] <- as.numeric(difftime(end.time, start.time, units='mins'))
    roc[["preda.best.case"]] <- roc.auc2(ordg$best_case_order, cancerGenes)
    roc[["preda.worst.case"]] <- roc.auc2(ordg$worst_case_order, cancerGenes)
    ordered.genes[["preda.best.case"]] <- ordg$best_case_order        
    ordered.genes[["preda.worst.case"]] <- ordg$worst_case_order            
  }
  
  # Run methods below only when we have two-group comparison setup
  if (length(unique(Labels)) == 2) {

    # DRI-SAM
    available.methods <- c(available.methods, "DRI.ss")    
    if (!is.null(methods) && ("DRI.ss" %in% methods) && evaluate) {
      message("DRI.ss")
      start.time <- Sys.time()      
      ordg <- test.geneorder.dri.sam(ge=ge2, cn=cn2, Labels, nperm=nperm, transform.type="standardize", version=version)
      end.time <- Sys.time()
      runtime[["DRI.ss"]] <- as.numeric(difftime(end.time, start.time, units='mins'))      
      roc[["DRI.ss"]] <- roc.auc2(ordg, cancerGenes)
      ordered.genes[["DRI.ss"]] <- ordg
    }

    available.methods <- c(available.methods, "DRI.srank")
    if (!is.null(methods) && ("DRI.srank" %in% methods) && evaluate) {
      message("DRI.srank")
      start.time <- Sys.time()      
      ordg <- test.geneorder.dri.sam(ge=ge2, cn=cn2, Labels, nperm=nperm, transform.type="rank", version=version)
      end.time <- Sys.time()
      runtime[["DRI.srank"]] <- as.numeric(difftime(end.time, start.time, units='mins'))      
      roc[["DRI.srank"]] <- roc.auc2(ordg, cancerGenes)
      ordered.genes[["DRI.srank"]] <- ordg      
    }
    
    available.methods <- c(available.methods, "DRI.sraw")    
    if (!is.null(methods) && ("DRI.sraw" %in% methods) && evaluate) {
      message("DRI.sraw")
      start.time <- Sys.time()      
      ordg <- test.geneorder.dri.sam(ge=ge2, cn=cn2, Labels, nperm=nperm, transform.type="raw", version=version)
      end.time <- Sys.time()
      runtime[["DRI.sraw"]] <- as.numeric(difftime(end.time, start.time, units='mins'))      
      roc[["DRI.sraw"]] <- roc.auc2(ordg, cancerGenes)
      ordered.genes[["DRI.sraw"]] <- ordg      
    }

  }

  return(list(auc = sapply(roc, function(x) {x$auc}), roc = roc, 
              runtime = runtime, ordered.genes = ordered.genes, 
	      cancerGenes = cancerGenes, available.methods = available.methods))
  
}
