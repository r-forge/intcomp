process.copynumber <- function (cn.raw, cn.seg, probespanCN, prior = "all", organism = "other") {
  
  require(CGHcall)
     
  if(length(cn.seg) > 0){ # segmented data given
   
    # If segmented data is given, just convert it into cghSeg format
    
    annots <- data.frame(name = rownames(cn.seg$data),
	                 Chromosome = as.integer(cn.seg$info$chr),
			 Start = as.integer(cn.seg$info$loc) - probespanCN,
                         End = as.integer(cn.seg$info$loc) + probespanCN)
    featureData_cn <- new("AnnotatedDataFrame", data = annots)
    rownames(featureData_cn@data) <- rownames(cn.seg$data)
    
    # if segmented data is given, we assume that the user has already
    # taken care of necessary preprocessing steps
    cgh.seg <- new('cghSeg', segmented = as.matrix(cn.seg$data), copynumber = as.matrix(cn.seg$data), featureData = featureData_cn)

  } else if(length(cn.seg) == 0){
  
    # Segment the data if segmented data is not available
    annots <- data.frame(name = rownames(cn.raw$data),
	                 Chromosome = as.integer(cn.raw$info$chr),
			 Start = as.integer(cn.raw$info$loc) - probespanCN,
                         End = as.integer(cn.raw$info$loc) + probespanCN)
			 
    cgh.raw <- make_cghRaw(cbind(annots, as.data.frame(cn.raw$data)))
    cgh.pre <- preprocess(cgh.raw, nchrom=length(unique(cgh.raw@featureData@data$Chromosome)))
    cgh.nor <- normalize(cgh.pre) # note: normalization median / smooth / none; using default: median
    cgh.seg <- segmentData(cgh.nor)
  }
      
  # Copy number calling
  print("cghcall")
  cgh.psn <- postsegnormalize(cgh.seg)
  cgh.cal <- CGHcall(cgh.psn, prior = prior, organism = organism) # number of copy number states can be changed with nclass parameter

  cgh <- ExpandCGHcall(cgh.cal, cgh.psn)
  
  return(cgh)
  
}
