# ge            list object contain a data matrix with gene expression data and a matrix with annotations for gene expression data
# cn            object contain a data matrix with copy number data and a matrix with annotations for copy number data
# Labels        labels to distinguish between patient (1) and reference (0) samples
# cancerGenes   character vector containing the labels of the cancer genes
# nperm         number of permutations to perform
# methods       which methods should be performed? see list.methods() for list of available methods

test.geneorder.pipeline <- function (ge, 
                                      cn.raw = NULL,
				     cn.seg = NULL,
				     cn.call = NULL,				    
				     cghCall = NULL,
				     cancerGenes, 
				     nperm = 1e2, 
				     methods = NULL, 
				     callprobs = NULL, 
				     evaluate = TRUE, 
				     cn.default = "segmented",
             references ="both"
				     ) {

  #########################################################

  # Pick segmented/called CN data from cghCall object, if not
  # separately given in the arguments

  if (!is.null(cghCall)) {    
    # here we must assume that matching of ge/cn probes has already
    # been carried out and the rows in cn.seg/cn.call/ge correspond
    if (is.null(cn.seg)) {
      dat <- assayDataElement(cghCall, 'segmented')
      cn.info <- ge$info
      if (nrow(dat) == nrow(cn.info)) {
        cn.seg <- list(data = dat, info = cn.info)
        rownames(cn.seg$data) <- rownames(ge$data)
      }
    }
    if (is.null(cn.call)) {
      dat <- assayDataElement(cghCall, 'calls')
      cn.info <- ge$info
      if (nrow(dat) == nrow(cn.info)) {     
        cn.call <- list(data = dat, info = cn.info)
        rownames(cn.call$data) <- rownames(ge$data)
      }
    }
  }

  ###############################################################################

  # Set the default copy number data (raw/segmented/called) to be used
  # in the comparisons, expect for methods with specific requirements.  
  if (cn.default == "raw") {
      # By default, use raw copy number (cn.raw) with all methods except 
      # intCNGEan (cn.seg) and CNAmet (cn.call)
      cn <- cn.raw
  } else if (cn.default == "segmented") {
      # By default, use segmented copy number (cn.seg) with all methods except 
      # CNAmet (cn.call)
      cn <- cn.seg      
  } else if (cn.default == "called") {
      # By default, use called copy number (cn.call) with all methods
      # CNAmet (cn.call)
      cn <- cn.call
  }  
  if (!nrow(cn$data) == nrow(ge$data)) {stop("CN/GE need to be matched prior to the analysis!")}
 
  #########################################################
  
  roc <- list()
  runtime <- list()
  ordered.genes <- list()
  available.methods <- c()

  available.methods <- c(available.methods, "edira")
  if (!is.null(methods) && ("edira" %in% methods) && evaluate) {
      message("edira")
      start.time <- Sys.time()      
      ordg <- test.geneorder.edira(ge, cn, references = references)
      end.time <- Sys.time()
      runtime[["edira"]] <- as.numeric(difftime(end.time, start.time, units='mins'))
      roc[["edira"]] <- roc.auc(ordg, cancerGenes)
      ordered.genes[["edira"]] <- ordg
  }

  available.methods <- c(available.methods, "DRI.cp")
  if (!is.null(methods) && ("DRI.cp" %in% methods) && evaluate) {
    message("DRI.cp")
    start.time <- Sys.time()
    ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="pearson")
    end.time <- Sys.time()
    runtime[["DRI.cp"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["DRI.cp"]] <- roc.auc(ordg, cancerGenes)
    ordered.genes[["DRI.cp"]] <- ordg
  }
  
  available.methods <- c(available.methods, "DRI.cs")  
  if (!is.null(methods) && ("DRI.cs" %in% methods) && evaluate) {
    message("DRI.cs")
    start.time <- Sys.time()    
    ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="spearman")
    end.time <- Sys.time()
    runtime[["DRI.cs"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["DRI.cs"]] <- roc.auc(ordg, cancerGenes)
    ordered.genes[["DRI.cs"]] <- ordg
  }
  
  available.methods <- c(available.methods, "DRI.ct")
  if (!is.null(methods) && ("DRI.ct" %in% methods) && evaluate) {
      message("DRI.ct")
      start.time <- Sys.time()      
      ordg <- test.geneorder.dri.cor(ge, cn, nperm=nperm, meth="ttest")
      end.time <- Sys.time()
      runtime[["DRI.ct"]] <- as.numeric(difftime(end.time, start.time, units='mins'))      
      roc[["DRI.ct"]] <- roc.auc(ordg, cancerGenes)
      ordered.genes[["DRI.ct"]] <- ordg      
  }
 
  available.methods <- c(available.methods, "SIM.full")
  if (!is.null(methods) && ("SIM.full" %in% methods) && evaluate) {
    message("SIM.full")
    start.time <- Sys.time()
    ordg <- test.geneorder.sim(ge, cn, meth = "full", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = seq(unique(ge$info$chr)))
    end.time <- Sys.time()
    runtime[["SIM.full"]] <- as.numeric(difftime(end.time, start.time, units='mins'))
    roc[["SIM.full"]] <- roc.auc(ordg, cancerGenes)
    ordered.genes[["SIM.full"]] <- ordg    
  }
  
  available.methods <- c(available.methods, "SIM.window")  
  if (!is.null(methods) && ("SIM.window" %in% methods) && evaluate) {
    message("SIM.window")
    start.time <- Sys.time()    
    ordg <- test.geneorder.sim(ge, cn, meth = "window", runname = paste("simtest-", abs(rnorm(1)), sep = ""), regs = seq(unique(ge$info$chr)), win = 1e6)
    end.time <- Sys.time()
    runtime[["SIM.window"]] <- as.numeric(difftime(end.time, start.time, units='mins'))
    roc[["SIM.window"]] <- roc.auc(ordg, cancerGenes)
    ordered.genes[["SIM.window"]] <- ordg
  }

  available.methods <- c(available.methods, "CNAmet")
  if (!is.null(methods) && ("CNAmet" %in% methods) && evaluate) {
     message("CNAmet")
     start.time <- Sys.time()     
     ordg <- test.geneorder.CNAmet(ge, cn = cn.call, nperm)
     end.time <- Sys.time()
     runtime[["CNAmet"]] <- as.numeric(difftime(end.time, start.time, units='mins'))     
     roc[["CNAmet"]] <- roc.auc(ordg, cancerGenes)
     ordered.genes[["CNAmet"]] <- ordg
  }
  
  available.methods <- c(available.methods, "intcngean")  
  if (!is.null(methods) && ("intcngean" %in% methods) && evaluate) {
    message("intCNGEan")
    start.time <- Sys.time()    
    if(is.null(callprobs)){
      ordg <- test.geneorder.intcngean(ge = ge,
                                       cghCall = cghCall,
                                       meth = "wmw",
                                       analysis.type = "univariate",
                                       nperm = nperm,
                                       pth = 0.1,
                                       match = TRUE)
    } else { #if(input == "simulations.equal.dimensions")
      ordg <- test.geneorder.intcngean(ge = ge, cghCall = cghCall, 
                                       meth = "wmw",
                                       analysis.type = "univariate",
                                       nperm = nperm, pth = 0.1, callprobs=callprobs, match=FALSE)
    }
    
    #if(input == "simulations.unequal.dimensions"){
    #  ordg <- test.geneorder.intcngean(ge=ge, cghCall=cghCall, 
    #                                   meth="wmw",
    #                                   analysis.type="univariate", 
    #                                   nperm = nperm, pth = 0.1, callprobs=NULL, match=TRUE)
    #
    #
    #}
    
    end.time <- Sys.time()
    runtime[["intCNGEan.wmw.univariate"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["intCNGEan.wmw.univariate"]] <- roc.auc(ordg, cancerGenes)
    ordered.genes[["intCNGEan.wmw.univariate"]] <- ordg
  }
  
  available.methods <- c(available.methods, "PMA.raw")  
  if (!is.null(methods) && ("PMA.raw" %in% methods) && evaluate) {
    message("PMA.raw")
    start.time <- Sys.time()    
    ordg <- test.geneorder.pma.rawscore(ge, cn)
    end.time <- Sys.time()
    runtime[["PMA.raw"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["PMA.raw"]]  <- roc.auc(ordg, cancerGenes)
    ordered.genes[["PMA.raw"]] <- ordg
  }
  
  available.methods <- c(available.methods, "pint")  
  if (!is.null(methods) && ("pint" %in% methods) && evaluate) {

    message("pint")
  
    start.time <- Sys.time()
    ordg <- test.geneorder.pint(ge, cn.raw, cn.seg)
    end.time <- Sys.time()
    runtime[["pint"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["pint"]] <- roc.auc(ordg, cancerGenes)
    ordered.genes[["pint"]] <- ordg    
  }

  available.methods <- c(available.methods, "OrtizEstevez")
  if (!is.null(methods) && ("OrtizEstevez" %in% methods) && evaluate) {
    message("OrtizEstevez")
    start.time <- Sys.time()    
    ordg <- test.geneorder.OrtizEstevez(ge, cn)
    end.time <- Sys.time()
    runtime[["OrtizEstevez"]] <- as.numeric(difftime(end.time, start.time, units='mins'))    
    roc[["OrtizEstevez"]] <- roc.auc(ordg, cancerGenes)
    ordered.genes[["OrtizEstevez"]] <- ordg        
  }                                      
          
  available.methods <- c(available.methods, "PREDA")	  
  if (!is.null(methods) && ("PREDA" %in% methods) && evaluate) {
    message("PREDA")
    start.time <- Sys.time()    
    ordg <- test.geneorder.preda(ge, cn, nperm=nperm, cancerGenes=cancerGenes,
        ge.qval.threshold=0.05, cn.qval.threshold=0.01, smoothMethod="spline",
        ge.smoothStatistic.threshold.up=0.5, ge.smoothStatistic.threshold.down=-0.5,
        cn.smoothStatistic.threshold.gain=0.1, cn.smoothStatistic.threshold.loss=-0.1, correction.method="fdr",
        chromosomes=unique(ge$info$chr))
    end.time <- Sys.time()    
    runtime[["preda"]] <- as.numeric(difftime(end.time, start.time, units='mins'))
    
    #ordg_num <- apply(ordg,2,as.numeric)
    #sort_ind <-function(x){
     # sort(x, index.return = TRUE)$ix
    #}
    
    #ordg_num_rank <- t(apply(ordg_num,1,sort_ind))
    #median_ranks <- apply(ordg_num_rank,2,median)
    
    roc_preda <- function(x){
       roc.auc(x, cancerGenes)$auc
    }

    roc_perms <- apply(ordg,1,roc_preda)
    choice <- which(roc_perms == median_int(roc_perms))
    roc[["preda"]] <- roc.auc(ordg[choice,], cancerGenes)
    ordered.genes[["preda"]] <- ordg[choice,]                
  }
  
  return(list(auc = sapply(roc, function(x) {x$auc}), roc = roc, 
  runtime = runtime, ordered.genes = ordered.genes, cancerGenes =
  cancerGenes, available.methods = available.methods))
  
}
