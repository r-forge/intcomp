test.geneorder.edira <- function(ge, cn, references){
  require(edira)
  # calculate edira
  edi <- edira_ratios(ge, cn, references)
  # order genes
  names(edi$test) <- rownames(ge$data)
  ordg <- names(sort(edi$test))  
  ordg
}
