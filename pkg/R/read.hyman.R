

read.hyman <- function (cdna, cgh, genenames, chrs = 1:22, xx) {

  # Load cdna, cgh, genenames 
  # preprocessed as in Berger et al., 2002
  # cdna: HymancdnaDataA.tab
  # cgh: HymancghDataA.tab
  # genenames: HymanAcc.mat
  # obtained from http://www.ece.ucsb.edu/pubs/ieee/index.shtml
  # Preprocessed as in Berger et al., 2002
  #cdna <- as.matrix(read.csv("/share/mi/data/pint10/Hyman2002/HymancdnaDataA.tab", sep = " ", header = FALSE))
  #cgh <- as.matrix(read.csv("/share/mi/data/pint10/Hyman2002/HymancghDataA.tab", sep = " ", header = FALSE))
  #genenames <- readLines("/share/mi/data/pint10/Hyman2002/HymanAcc.mat")
  
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

  # match the probes
  require(pint)
  tmp <- pint.match(geneExp, geneCopyNum, chrs = chrs)
       
  list(ge = tmp$X, cn = tmp$Y, Labels = NULL)

}
	    
	   
