roc.auc <- function (ordered.results, P) {

  # Calculate ROC/AUC value for ordered list, given the list P of true
  # positives
  
  # ordered results: best to worst
  # P: known positives

  # Compute area under curve
  rates <- roc(ordered.results, P)
  # integration: Compute step intervals and compute weighted sum of true positive rates in each interval.
  # note that interval length can be 0 if fpr does not change
  as.numeric(t(rates$fpr[-1]-rates$fpr[-length(rates$fpr)])%*%rates$tpr[-length(rates$fpr)])
  # Integration
  #sapply(names(rates$roc), function (nam){mean(na.omit(rates$roc$tpr[1:n]))})


}


roc.auc2 <- function (ordered.results, P) {

  # Calculate ROC/AUC value for ordered list, given the list P of true
  # positives
  
  # ordered results: best to worst
  # P: known positives

  # Compute area under curve
  rates <- roc(ordered.results, P)
  # integration: Compute step intervals and compute weighted sum of true positive rates in each interval.
  # note that interval length can be 0 if fpr does not change
  list(auc = as.numeric(t(rates$fpr[-1]-rates$fpr[-length(rates$fpr)])%*%rates$tpr[-length(rates$fpr)]), tpr = rates$tpr, fpr = rates$fpr)
}

roc <- function (ordered.results, P) {

  # Calculate ROC curve
    
  #ordered results: best to worst
  #P: known positives
  #output: true positive rate and false positive rate^M
		        
  #Check that all known positives are included in the original analysis i.e. ordered results list^M
  #if (!all(P %in% ordered.results)) {print("Warning: not all known positives are in the results list. Only included positives are used.")}
  positives <- P[P %in% ordered.results]
  #Number of retrieved known cytobands
  N<-length(ordered.results) #total number of samples
  Np<-length(positives) #number of positives
  Nn<-N-Np #number of negatives
  TP<-cumsum(ordered.results %in% positives)
  FP<-cumsum(!(ordered.results %in% positives))
  tpr <- TP/Np #TP/(TP + FN) = TP.simCCA / P
  fpr <- FP/Nn #FP/(FP + TN) = FP.simCCA / N.simCCA
  
  list(tpr=tpr,fpr=fpr)
}

	
