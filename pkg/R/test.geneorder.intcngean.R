test.geneorder.intcngean <- function (ge, cn=NULL, cn.raw, Labels=NULL, meth, analysis.type, probespanCN, probespanGE, nperm, pth, callprobs=NULL, prior="all", organism = "other", match="FALSE") {

  require(intCNGEan)
  if(length(cn) > 0){X.seg <- cn}
  X.raw <- cn.raw
  Y <- ge
  
  #print("--- CGH data ---")
  
  chrs <- as.character(X.raw$info[, "chr"])
  chrs[chrs == "23"] = "X"
  chrs[chrs == "24"] = "Y"  
  yprobes <- (chrs %in% c("X", "Y"))
  if (sum(yprobes) > 0) {
    X.raw$info <- X.raw$info[!yprobes,]
    X.raw$data <- X.raw$data[!yprobes,]
  }

  chrs <- as.character(Y$info[, "chr"])
  chrs[chrs == "23"] = "X"
  chrs[chrs == "24"] = "Y"  
  yprobes <- (chrs %in% c("X", "Y"))
  if (sum(yprobes) > 0) {
    warning("Removed probes in X and Y chromosomes, ", sum(yprobes), " probes in total out of ", 
    length(yprobes), " probes in the data.")
    Y$info <- Y$info[!yprobes,]
    Y$data <- Y$data[!yprobes,]
    if (length(cn) > 0) {
      X.seg$info <- X.seg$info[!yprobes,]
      X.seg$data <- X.seg$data[!yprobes,]
    }
  }
      
  ##################################################

  if(length(cn) > 0){ # segmented data given

    annots <- data.frame(name = rownames(X.seg$data),
                         Chromosome = as.integer(X.seg$info$chr),
                         Start = as.integer(X.seg$info$loc) - probespanCN,
                         End = as.integer(X.seg$info$loc) + probespanCN)
    featureData_cn <- new("AnnotatedDataFrame", data = annots)
    rownames(featureData_cn@data) <- rownames(X.seg$data)
  
    # if segmented data is given, we assume that the user has already
    # taken care of necessary preprocessing steps

    cghdata.seg <- new('cghSeg', segmented = as.matrix(X.seg$data), copynumber = as.matrix(X.seg$data), featureData = featureData_cn)

  } else if(length(cn) == 0){

    annots <- data.frame(name = rownames(X.raw$data),
                         Chromosome = as.integer(X.raw$info$chr),
                         Start = as.integer(X.raw$info$loc) - probespanCN,
                         End = as.integer(X.raw$info$loc) + probespanCN)

    cgh.tmp <- cghRaw(cbind(annots, as.data.frame(X.raw$data)))   
    cghdata.raw <- preprocess(cgh.tmp, nchrom = length(unique(cgh.tmp@featureData@data[, "Chromosome"])))  
    cghdata.seg <- segmentData(cghdata.raw)
  }

  #CghSeg <- new("cghSeg",featureData = featureData_cn, copynumber = as.matrix(copynumber(cghdata.raw)), segmented = as.matrix(copynumber(cghdata.seg)))
  print("cghcall")
  cghcall <- CGHcall(cghdata.seg, prior=prior, organism=organism)
  #cghcall <- CGHcall(cghdata.seg) # calling, name was cghSet, think prior
  print("callprobs")
  if(length(callprobs) > 0){
    probgain(cghcall) <- callprobs$gain
    probloss(cghcall) <- callprobs$loss
    probnorm(cghcall) <- callprobs$norm
  }
 
  #--- MRNA data ---
  # Add start and end positions
  ychr <- Y$info$chr
  ychr <- gsub("X",23,ychr)
  ychr <- gsub("Y",24,ychr)

  #annots <- data.frame(name = as.character(rownames(Y$data)),
  #                 Chromosome = as.numeric(as.character(ychr)),
  #                 Start = as.integer(Y$info$loc) - probespanGE,
  #                 End =  as.integer(Y$info$loc) + probespanGE)
  #featureData_ge <- new("AnnotatedDataFrame", data = annots)
  #rownames(featureData_ge@data) <- rownames(Y$data)

  #mrnaSet <- new("ExpressionSet", exprs = as.matrix(Y$data), featureData = featureData_ge)
  ############################
  print("MRNA")
  df <- data.frame(Chromosome = as.numeric(as.character(ychr)),
                   Start = as.integer(Y$info$loc) - probespanGE,
                   End = as.integer(Y$info$loc) + probespanGE)

  featureData_ge <- new("AnnotatedDataFrame",
                         varMetadata = data.frame(labelDescription = c("Chromosomal position","Basepair position start","Basepair position end"),
                         row.names = c("Chromosome", "Start", "End")),dimLabels = c("featureNames", "featureColumns"), data = df)

  featureNames(featureData_ge) <- rownames(Y$data)
  mrnaSet <- new("ExpressionSet", exprs = as.matrix(Y$data), featureData = featureData_ge)

  # Match CN and GE
  if(match == TRUE){
    print(" Match data sets")
    matched <- intCNGEan.match(cghcall, mrnaSet, GEbpend = "yes", CNbpend = "yes")
    cghcall <- matched$CNdata.matched
    mrnaSet <- matched$GEdata.matched
  }
  # already matched, so just read from the data
  #matched <- list(CNdata.matched = cghcall, GEdata.matched = mrnaSet) # was uncommented with Leo
 
  #meth = "wmw"    # weighted Mann-Whitney
  #meth = "wcvm"   # weighted Cramer-Von Mises; failed

  print(" Tune model parameters")

  #tuned  <- intCNGEan.tune(cghcall, mrnaSet, test.statistic=meth, ngenetune = 250, nperm_tuning = 250, minCallProbMass = 0.01)
  #tuned  <- intCNGEan.tune(matched$CNdata.matched, matched$GEdata.matched, test.statistic = meth)
  #tuned  <- intCNGEan.tune(cghcall, mrnaSet, test.statistic = meth) # try this if does not work with parameters
  tuned  <- intCNGEan.tune(cghcall, mrnaSet, test.statistic = meth, ngenetune = 250, nperm_tuning = 250, minCallProbMass = 0.01)
 
  # Calculate the models
  #for (analysis.type in c("univariate")) { # regional caused error

  tested <- intCNGEan.test.mod(tuned, analysis.type = analysis.type, test.statistic = meth, nperm = nperm, eff.p.val.thres = pth)

  # If this works, then try mod version
  #tested <- intCNGEan.test(tuned, analysis.type = analysis.type, test.statistic = meth, nperm = nperm, eff.p.val.thres = pth)

  #ind_pval <- which(colnames(tested)=="raw.p")
  #ind_names <- which(colnames(tested)=="name")
  #pvals <- unlist(tested[,ind_pval])
  #names <- as.character(unlist(tested[,ind_names]))
  ## Order the genes
  #ordg <- names[sort(pvals, index.return= TRUE)$ix]
  #ordg  

  print(" Order the genes")
  ordg <- rownames(tested)[order(tested[["raw.p"]], decreasing = FALSE)]
  ordg
}

  
