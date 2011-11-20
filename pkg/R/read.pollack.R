# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.


read.pollack <- function (chrs = 1:22, dat, clone2geneid) {

  #http://www.pnas.org/content/suppl/2002/09/23/162471999.DC1/4719CopyNoGeneDatsetLegend.html
		       
  cn <- dat[, 3:48]
  colnames(cn) <- unname(sapply(colnames(cn), function(s) { x <- paste(unlist(strsplit(s, "\\."))); paste(x[1:min(length(x),2)], collapse = "")}))
  rownames(cn) <- as.character(dat[,1])
    
  ge <- dat[, 49:89]
  colnames(ge) <- unname(unlist(sapply(colnames(ge), function(s) {paste(unlist(strsplit(s, "\\."))[1:2], collapse = "")})))
  rownames(ge) <- as.character(dat[,1])

  # Impute missing values
  cn <- as.matrix(impute(cn))
  ge <- as.matrix(impute(ge))

  # Normalize (unit variance and zero mean for each column)
  #cn <- centerData(unitscale(cn))
  #ge <- centerData(unitscale(ge))
      
  # Convert Pollack clone IDs to GeneIDs
  # Clone names in Pollack data (IMAGE:*****) refer to known
  # clones. Mapped IMAGE:**** symbols to EntrezGeneIDs with the tool
  # at http://smd.stanford.edu/cgi-bin/source/sourceBatchSearch
  # resulting in 'PollackClone2GeneID.tab' file (stored in data
  # directory in data(pollack) clone2geneid object). We omitted genes that are not in EntrezGene/UniGene,
  # and selected option 'Show all IDs if in multiple Genes/Clusters'.

  # Pick clones that have corresponding GeneID       
  tab <- clone2geneid                                      
  clones <- as.character(dat[,1])  
  # Remove clones that match with multiple genomic locations
  tab <- tab[-grep("In multiple ClusterIDs",tab[,2]),]
    
  # Remove empty geneIDs
  tab <- tab[!tab["GeneID"] == "",]
  gids <- as.character(tab[, 2])
  names(gids) <- as.character(tab[, 1])
	    
  # Return clone-to-geneid mapping, remove NAs and duplicates
  clone2gid <- gids[clones]
  clone2gid <- clone2gid[!is.na(clone2gid)]
  clone2gid <- clone2gid[!duplicated(clone2gid)]
		    

  # Pick clones for which geneID was found
  ge <- ge[names(clone2gid),]
  cn <- cn[names(clone2gid),]
  # Name clones by GeneID
  rownames(ge) <- rownames(cn) <- clone2gid

  # Get chromosomal locations for the GeneIDs
  info <- get.entrez.info(rownames(ge))

  # take those genes (clones) for which location found
  genes <- intersect(rownames(ge), rownames(info))  

  # Pick common samples and in the same order and include only
  # probes with location information
  commons <- intersect(colnames(ge), colnames(cn))
  cn <- cn[genes, commons]
  ge <- ge[genes, commons]
    
  geneExp <- list(data = ge, info = info) 
  geneCopyNum <- list(data = cn, info = info)

  # Note we can match ge/cn already before segmentation
  # as they have the same resolution
  require(pint)                                             
  tmp <- pint.match(geneExp, geneCopyNum, chrs = chrs)

  # Segmentation/calling for copy number data
  CN.raw <- tmp$Y
  cgh <- process.copynumber(CN.raw)
  #CN.seg <- list(data = assayDataElement(cgh, 'segmented'), info = CN.raw$info)
  #CN.call <- list(data = assayDataElement(cgh, 'calls'), info = CN.raw$info)
  #rownames(CN.call$data) <- rownames(CN.seg$data) <- rownames(geneExp$data)

  list(ge = tmp$X, cn.raw = tmp$Y, cghCall = cgh)
}
