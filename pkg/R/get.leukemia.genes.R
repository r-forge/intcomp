
get.leukemia.genes <- function (stratton.cancer.gene.census) {

  tab <- stratton.cancer.gene.census

  g <- NULL
  g <- c(g, get.stratton.genes(tumor.type = "leukemia", tab))
  g <- c(g, get.stratton.genes(tumor.type = "ALL", tab))
  g <- c(g, get.stratton.genes(tumor.type = "AML", tab))
  g <- c(g, get.stratton.genes(tumor.type = "CLL", tab))
  g <- c(g, get.stratton.genes(tumor.type = "APL", tab))
  g <- c(g, get.stratton.genes(tumor.type = "NHL", tab))
  g <- c(g, get.stratton.genes(tumor.type = "CMML", tab))
  g <- c(g, get.stratton.genes(tumor.type = "ALCL", tab))
  g <- c(g, get.stratton.genes(tumor.type = "DLBCL", tab))
  g <- c(g, get.stratton.genes(tumor.type = "MDS", tab))
  g <- c(g, get.stratton.genes(tumor.type = "AEL", tab))
  g <- c(g, get.stratton.genes.exactmatch(tumor.type = "AL", tab))
  g <- unique(na.omit(g))

  g
}



get.stratton.genes <- function (tumor.type = NULL, tab) {

  # List known cancer genes in Stratton's cancer gene list
  #tab <- read.csv("/share/mi/data/cancer_genes/Table_1_full_2010_02_11.txt", sep = "\t")
  # data(stratton.cancer.gene.census); tab <- stratton.cancer.gene.census
  genes <- as.character(tab[, "GeneID"])

  if (!is.null(tumor.type)) {
    somatic <- which(sapply(sapply(tab[,8], function(x){grep(tumor.type,as.character(x))}), length) > 0)
    germ    <- which(sapply(sapply(tab[,9], function(x){grep(tumor.type,as.character(x))}), length) > 0)
    inds    <- union(somatic, germ)
    genes   <- genes[inds]
  }
  
  unique(genes)

}

get.stratton.genes.exactmatch <- function (tumor.type = NULL, tab) {

  # Here find only exact matches of the tumor.type string
  
  # List known cancer genes in Stratton's cancer gene list
  # Mike Stratton's list (From JHO 3/2010)
  # data(stratton.cancer.gene.census); tab <- stratton.cancer.gene.census  

  genes <- as.character(tab[, "GeneID"])

  if (!is.null(tumor.type)) {
    somatic <- which(sapply(sapply(tab[,8], function(x){which(tumor.type == as.character(x))}), length) > 0)
    germ <- which(sapply(sapply(tab[,9], function(x){which(tumor.type == as.character(x))}), length) > 0)
    inds <- union(somatic, germ)
    genes <- genes[inds]
  }

  unique(genes)

}


