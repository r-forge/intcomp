test.geneorder.edira <- function(ge, cn, Labels=NULL){
  require(edira)
  # calculate edira
  edi <- edira_ratios(cbind(ge$data,ge$info$chr,(ge$info$end+ge$info$start)/2), cbind(cn$data,cn$info$chr,(cn$info$end+cn$info$start)/2))
  # order genes
  names(edi$test) <- rownames(ge$data)
  ordg <- names(sort(edi$test))  
  ordg
}
