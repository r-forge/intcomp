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
}
