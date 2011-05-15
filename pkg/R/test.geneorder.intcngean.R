test.geneorder.intcngean <- function (ge, cn, Labels=NULL, meth, analysis.type, nperm, pth, callprobs=NULL) {

  # NOTE: assumes cghCall object in input!
  mrnaSet <- process.ge(ge, probespanGE = 16) # Convert ge into expressionSet object 

  require(intCNGEan)
	  
  print("callprobs")
  if(length(callprobs) > 0){
    probgain(cn) <- callprobs$gain
    probloss(cn) <- callprobs$loss
    probnorm(cn) <- callprobs$norm
  }
 
  # meth = "wmw"    # weighted Mann-Whitney
  # meth = "wcvm"   # weighted Cramer-Von Mises; failed
  # analysis.type regional caused error

  print(" Tune model parameters")
  tuned  <- intCNGEan.tune(cn, mrnaSet, test.statistic = meth, ngenetune = 250, nperm_tuning = 250, minCallProbMass = 0.01)
  #tuned  <- intCNGEan.tune(cn, mrnaSet, test.statistic=meth, ngenetune = 250, nperm_tuning = 250, minCallProbMass = 0.01)
  #tuned  <- intCNGEan.tune(matched$CNdata.matched, matched$GEdata.matched, test.statistic = meth)
  #tuned  <- intCNGEan.tune(cn, mrnaSet, test.statistic = meth) # try this if does not work with parameters
   
   # Calculate the models
  tested <- intCNGEan.test.mod(tuned, analysis.type = analysis.type, test.statistic = meth, nperm = nperm, eff.p.val.thres = pth)

  print(" Order the genes")
  ordg <- rownames(tested)[order(tested[["raw.p"]], decreasing = FALSE)]
  ordg
}

  
