list.methods <- function (Labels = c(1,2)) {
  res <- test.geneorder.pipeline(cn.raw = NULL, cn.seg = NULL, cn.call = NULL, 
          Labels = Labels, cancerGenes = NULL, evaluate = FALSE)
  res$available.methods
}