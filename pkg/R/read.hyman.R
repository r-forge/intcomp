# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.


read.hyman <- function (cdna, cgh, genenames, chrs = 1:22, xx, useSegmentedData = FALSE, remove.duplicates = TRUE) {

 #for (f in list.files("~/local/Rpackages/intcomp/intcomp/intcomp/pkg/R/", full.names = TRUE, pattern = ".R$")) {source(f)}
 #xx = as.list(org.Hs.egALIAS2EG); useSegmentedData = TRUE; remove.duplicates = FALSE; chrs <- 1:22

  # Load cdna, cgh, genenames 
  # preprocessed as in Berger et al., 2002
  # cdna: HymancdnaDataA.tab
  # cgh: HymancghDataA.tab
  # genenames: HymanAcc.mat
  # obtained from http://www.ece.ucsb.edu/pubs/ieee/index.shtml
  # Preprocessed as in Berger et al., 2002
  # cdna <- as.matrix(read.csv("HymancdnaDataA.tab", sep = " ", header = FALSE))
  # cgh <- as.matrix(read.csv("HymancghDataA.tab", sep = " ", header = FALSE))
  # genenames <- readLines("HymanAcc.mat")
  
  # Remove rows with no gene name
  keep <- !genenames == ""   
  cdna <- cdna[keep,]
  cgh <- cgh[keep,]
  rownames(cdna) <- rownames(cgh) <- genenames[keep]	      
  genes.uniq <- unique(genenames[keep])
			
  # Pick only those geneids where a unique symbol-geneid match is found
  gids <- sym2gid(genes.uniq, xx, uniq = FALSE) # convert to gids
  gids <- sapply(gids[sapply(gids, function(x){length(x[[1]])}) == 1], function(x){x[[1]]})

  # Now pick the rows of the data matrices that correspond to a unique GeneID
  genes <- intersect(rownames(cdna), names(gids))  
  cdna <- cdna[genes,]
  cgh <- cgh[genes,]
      
  #name rows by GeneIDs instead of symbols
  rownames(cdna) <-  rownames(cgh) <- unname(gids[rownames(cdna)])
	  
  # Get chromosomal locations for GeneIDs
  info <- get.entrez.info(rownames(cdna))
	        		
  # take those genes for which location found
  genes <- intersect(rownames(cdna), rownames(info))  
  cdna <- cdna[genes,]
  cgh <- cgh[genes,]
  info <- info[genes,]
	  
  # Organize
  geneExp <- list(data = cdna, info = info) 
  geneCopyNum <- list(data = cgh, info = info)

  # match the probes using tools in pint package
  # note: cn and ge already at the same resolution,
  # no need to segment cn before matching!
  require(pint)
  tmp <- pint.match(geneExp, geneCopyNum, chrs = chrs, useSegmentedData = useSegmentedData, remove.duplicates = remove.duplicates)

  # Segmentation/calling for copy number data
  CN.raw <- tmp$Y
  if (any(duplicated(rownames(CN.raw$data)))) {
    nams <- apply(cbind(rownames(CN.raw$data), 1:nrow(CN.raw$data)), 1, function (x) {paste(x, collapse = "-rowid")})
    rownames(CN.raw$data) <- nams
    rownames(CN.raw$info) <- nams
  } 
  cgh <- process.copynumber(CN.raw)

  #list(ge = tmp$X, cn = tmp$Y, Labels = NULL)
  list(ge = tmp$X, cn.raw = tmp$Y, cghCall = cgh)

}
	    
	   
