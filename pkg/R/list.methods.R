list.methods <- function (Labels = c(1,2)) {
  res <- test.geneorder.pipeline(ge = list(data = matrix(rnorm(9),3,3)),
                                 cn.raw = list(data = matrix(rnorm(9),3,3)), 
          cancerGenes = NULL, evaluate = FALSE, cn.default = "raw")
  res$available.methods
}