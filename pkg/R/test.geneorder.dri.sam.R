test.geneorder.dri.sam <- function (ge, cn, Labels=NULL, nperm = 1e2, transform.type, version) {
  require(DRI)
  
  DNA.data <- as.matrix(cn$data[,Labels==1])
  RNA.data <- as.matrix(ge$data[,Labels==1])
  DNA.norm.data <- as.matrix(cn$data[,Labels==0])
  RNA.norm.data <- as.matrix(ge$data[,Labels==0])
  Labels.sam <- c(rep(1,ncol(DNA.data)),rep(2,ncol(DNA.norm.data)))
  if(version=="normal"){
  auswahl <- sample(1:nrow(DNA.data), nperm)
    dr_ttest_results <- function(z){
        permute_dr_sam <- function(pos_ge){
           drsam(DNA.data=matrix(rep(c(DNA.data[z,],DNA.norm.data[z,]),2), nrow=2, byrow=TRUE), RNA.data=matrix(rep(c(RNA.data[pos_ge,],RNA.norm.data[pos_ge,]),2), nrow=2, byrow=TRUE), labels=Labels.sam, transform.type = transform.type)$test.summed[1]
        }
    obs <- permute_dr_sam(z)             
    verteilung_dr_sam <- as.numeric(sapply(auswahl,permute_dr_sam))
    pwert_dr_sam <- length(which(sort(abs(verteilung_dr_sam)) >= abs(obs)))/length(verteilung_dr_sam)               
    pwert_dr_sam
    }
  results <- sapply(1:nrow(DNA.data),dr_ttest_results)  
  }
  
  if(version=="approx"){
  permutation <- sample(1:nrow(DNA.data), nrow(DNA.data))
  obs <- drsam(DNA.data, RNA.data, Labels, transform.type = transform.type)
  distribution_dr_sam <- drsam(DNA.data, RNA.data[permutation,], labels=Labels.sam, transform.type = transform.type)
        calculate_pvalue <- function(z){
            length(which(sort(distribution_dr_sam) >= obs[z]))/length(distribution_dr_sam)                       
        }
  results <- sapply(1:length(obs),calculate_pvalue)  
  }
  
  # Gene names ordered by empirical p-value
  names(results) <- rownames(DNA.data)
  ordg <- names(results)[order(results)]
  ordg 
}
  
