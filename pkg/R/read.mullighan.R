read.mullighan <- function(chrs = 1:22, location.table, cnAnnotation, cnData, cn, eSet) {
  
  require(affy)

  ### Expression data ###  
  # The probeset IDs are "ENTREZID_at"
  #load("/share/mi/data/CMG/integration-review-2010/Mullighan/eSet.RData")
  #load(paste(data.path, "eSet.RData", sep = "/"))  
  X <- exprs(eSet)
  medExprs <- apply(exprs(eSet), 1, median) # median expression
  X <- exprs(eSet) <- exprs(eSet) - medExprs
  
  # Annotations for the genes
  #data.path + hgu133ahsentrezgcdf_13.0.0.tar.gz
  #install.packages("hgu133ahsentrezg.db_13.0.0.tar.gz")
  #location information for TSSs
  #from annotation package "hgu133ahsentrezg.db_13.0.0.tar.gz"

  sets <- removeAFFX(rownames(X))
  xinfo <- get.TSS.locations.for.entrezids(sets, location.table)
  xdata <- X[rownames(xinfo),]

  #?hgu133ahsentrezgCHRLOC
  #Mappings were based on data provided by: UCSC Genome
  #Bioinformatics (Homo sapiens)
  #ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19 With a date stamp
  #from the source of: 2009-Jul5

  #?hgu133ahsentrezgENTREZID
  #Mappings were based on data provided by: Entrez Gene
  #ftp://ftp.ncbi.nlm.nih.gov/gene/DATA With a date stamp from the
  #source of: 2010-Mar1

  #load(paste(data.path, "rawCn.RData", sep = "/"))  
  #Y <- NULL
  #for (i in chrs) {Y <- rbind(Y, cn[[i]])}

  ### Expression data ###
  
  # Segmented CN data: mcr x patients matrix
  # "cnAnnotation" "cnData"
  #load(paste(data.path, "cnMCRs.RData", sep = "/"))  
  # Each array segmented separately, no MCRs
  #load(paste(data.path, "segmentedCn.RData", sep = "/"))

  # Raw CN data; required by intCNGEan
  #load(paste(data.path, "rawCn.RData", sep = "/")) # cn
  #cn.raw <- as.matrix(ldply(cn, function (x) {x})[,-1])
  cn.raw <- NULL
  info.raw <- NULL
  for (cni in names(cn)) {
    print(cni)
    cn.raw <- rbind(cn.raw, cn[[cni]])
    chr <- as.numeric(cni) # list element names are chromosome names
    locs <- as.numeric(rownames(cn[[cni]])) # feature names are chromosome names
    info.raw <- rbind(info.raw, cbind(rep(chr, length(locs)), locs))
    cn[[cni]] <- NULL    # free some memory

  }
  colnames(info.raw) <- c("chr", "loc")
  rm(cn)
  rownames(info.raw) <- rownames(cn.raw) <- as.character(seq(1:nrow(cn.raw)))
  
  # format column names so that they are comparable
  # i.e. remove .CEL ending from ge data columns
  colnames(xdata) <- unname(sapply(colnames(xdata), function(x) {strsplit(x, "\\.")[[1]][[1]]}))
  
  # Form data objects
  ge <- list(data = as.matrix(xdata), info = as.data.frame(xinfo))
  cn <- list(data = as.matrix(cnData), info = as.data.frame(cnAnnotation))
  cn.raw <- list(data = as.matrix(cn.raw), info = as.data.frame(info.raw))
  
  # Match rows
  library(pint)
  matched <- pint.match(ge, cn, max.dist = 1e7, chrs = chrs)
  
  list(ge = matched$X, cn = matched$Y, labels = unname(sapply(colnames(matched$X$data), function(x){strsplit(x,"_")[[1]][[1]]})), cn.raw = cn.raw)

}


get.TSS.locations.for.entrezids <- function (sets, location.table) {

  chrlist <- location.table[["chr"]]
  geneStart <- location.table[["geneStart"]]
  geneEnd <- location.table[["geneEnd"]]

  sets <- intersect(intersect(names(which(!is.na(geneStart))), names(which(!is.na(geneEnd)))), names(which(!is.na(chrlist))))
  chrlist <- chrlist[sets]
  geneStart <- geneStart[sets]
  geneEnd <- geneEnd[sets]
  
   tss.list <- c()
  multi.tss <- c()

  for (i in 1:length(sets)) {
      chr <- chrlist[[i]]
    start <- geneStart[[i]]
      end <- geneEnd[[i]]
    
    if (length(chr) > 1) {stop(paste("Multiple chrs mapped! Set ", sets[[i]]))}
    # If all positions negative, then gene is on neg. strand
    # and end position gives TSS
    if (all(c(sign(start) < 0, sign(end) < 0))) {
      tss <- unique(end)        
    } else if (all(c(sign(start) > 0, sign(end) > 0))) {
    # If all positions positive, then gene is on pos. strand
      tss <- unique(start)        
    } else {
      #stop("Both neg. and pos. location info given?")
      tss <- Inf
    }
    # If multiple TSSs take mean
    if (length(tss) > 1) {
      tss <- mean(tss)
      multi.tss <- c(multi.tss, i)
    }
    tss.list[[i]] <- abs(tss)
  }   
  message(paste("Proportion of sets with multiple TSSs: ", length(multi.tss)/length(sets)))

  # Sanity check
  if (!length(tss.list) == length(sets)) {stop("Error 1")}
  if (!length(chrlist) == length(sets)) {stop("Error 2")}
  
  # Remove NAs (gene readings on both strands; therefore no unique TSS info)
  inds <- !is.infinite(tss.list)
  sets <- sets[inds]
  tss.list <- tss.list[inds]
  chrlist <- chrlist[inds]

  message(paste("Proportion of sets with readings on both strands: ", sum(!inds)/length(sets)))  

  X.info <- cbind(chrlist, tss.list)
  colnames(X.info) <- c("chr","loc")
  rownames(X.info) <- sets
  
  return(X.info)
  
}

