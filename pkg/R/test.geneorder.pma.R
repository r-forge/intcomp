

test.geneorder.pma <- function (ge, cn, Labels = NULL, nperm = 1e2) {

  # Note that this method does not utilize Label information, option added for compatibility

  # Calculate empirical pvalues with permutations to order the genes

  ## use type="ordered"
  ## since both gene expression and copy number data sets are
  ## ordered by their chromosomal location

  require(PMA)
  
  tab <- get.pma.table(cn, ge, typex = "ordered", typez = "ordered")
        
  #############################################
	   
  # Calculate scores for random permutations and form a table of scores
  rand.scores <- array(NA, dim = c(nrow(ge$data), nperm))
  rownames(rand.scores) <- rownames(ge$data)
  for (i in 1:nperm) {
    # Random permutation for expression probes
    ge.rand <- ge
    ge.rand$data <- ge$data[sample(nrow(ge$data)),]
    rtab <- get.pma.table(cn, ge.rand, typex = "ordered", typez = "ordered")
      rand.scores[, i] <- rtab[, "score"]
    }

  # Empirical pvalues: high score is better, so calculate proportion
  # of values that are worse than random scores to get pvalues
  pvals <- apply((tab[, "score"] - rand.scores) < 0, 1, mean)
	    
  ############################################
	       
  #tab.ordered <- esort(tab, -score)
  tab.ordered <- tab[order(pvals), ]
		     
  # order genes
  ordg <- tab.ordered[, "probeid"]

  ordg
}


