test.geneorder.dri.cor <- function (ge, cn, Labels = NULL, nperm, meth, version) {
  require(DRI)
  
  DNA.data <- as.matrix(cn$data)
  RNA.data <- as.matrix(ge$data)  
  
  # DR-Correlate analysis to find genes with correlated DNA/RNA measurements
  if(version=="normal"){
  auswahl <- sample(1:nrow(DNA.data), nperm, replace = TRUE)
    if(meth=="pearson" | meth=="spearman"){
  # Get null distribution based on gene permutations
        dr_corr_results <- function(z){
            permute_dr_corr <- function(pos_ge){
                drcorrelate( matrix(DNA.data[z,], nrow=1), matrix(RNA.data[pos_ge,], nrow=1), method = meth)
            }
        obs <- permute_dr_corr(z)
        distribution_dr_corr <- as.numeric(sapply(auswahl,permute_dr_corr ))
        pvalue_dr_corr <- length(which(sort(distribution_dr_corr) >= obs))/length(distribution_dr_corr)                       
        pvalue_dr_corr
        }
    results <- sapply(1:nrow(DNA.data),dr_corr_results) 
    }
    if(meth=="ttest"){
    # Get null distribution based on gene permutations
        dr_ttest_results <- function(z){
            permute_dr_ttest <- function(pos_ge){
                drcorrelate(DNA=matrix(rep(DNA.data[z,],2), nrow=2, byrow=TRUE), RNA=matrix(rep(RNA.data[pos_ge,],2), nrow=2, byrow=TRUE), method = "ttest")[1]
            }
        obs <- permute_dr_ttest(z)               
        distribution_dr_ttest <- as.numeric(sapply(auswahl,permute_dr_ttest))
        pvalue_dr_ttest <- length(which(sort(distribution_dr_ttest) >= obs))/length(distribution_dr_ttest)                       
        pvalue_dr_ttest
        }
    results <- sapply(1:nrow(DNA.data),dr_ttest_results)  
    }
  }
  
  if(version=="approx"){
  permutation <- sample(1:nrow(DNA.data), nrow(DNA.data))
    if(meth=="pearson" | meth=="spearman"){
    # Get null distribution based on sample permutations
        obs <- drcorrelate(DNA=DNA.data, RNA=RNA.data, method = meth)
        distribution_dr_corr <- drcorrelate(DNA.data, RNA.data[permutation,], method = meth)
        calculate_pvalue <- function(z){
            length(which(sort(distribution_dr_corr) >= obs[z]))/length(distribution_dr_corr)                       
        }
    results <- sapply(1:length(obs),calculate_pvalue) 
    }
    if(meth=="ttest"){
    # Get null distribution based on sample permutations
        obs <- drcorrelate(DNA=DNA.data, RNA=RNA.data, method = "ttest")[1]
        distribution_dr_ttest <- drcorrelate(DNA=DNA.data, RNA=RNA.data[permutation,], method = "ttest")
        calculate_pvalue <- function(z){
            length(which(sort(distribution_dr_ttest) >= obs[z]))/length(distribution_dr_ttest)                       
        } 
    results <- sapply(1:nrow(DNA.data),calculate_pvalue)  
    }
  }
  # Gene names ordered by empirical p-value
  names(results) <- rownames(DNA.data)
  ordg <- names(results)[order(results)]
  ordg 
}
