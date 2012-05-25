# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.


#### Generates four data sets, one cases and one control data set for each copy number and gene expression (no ratios of cases and controls),
#### based on the input of one real copy number data set and one real gene expression data set, both supposedly from cases.
#### The simulation setup roughly follows the one of Schäfer et al. (2009): Integrated analysis of copy number alterations and gene expression: a bivariate assessment of equally directed abnormalities.
#### A joint distribution for the CN and GE values is assumed.

# GE_data      a matrix of gene expression data, without annotation
# CN_data      a matrix of copy number data, without annotaion
# Inner        indices of inner grid points for simulation of the copy number data
# Outer        indices of outer grid points for simulation of the copy number data
#              (one mixture component of the normal mixture distribution has as copy number coordinate of its mean an inner grid point, the other mixture component has as copy number coordinate of its mean an outer grid point)
# probs_GE     the quantiles used to calculate the GE grid points
# probs_CN     the quantiles used to calculate the CN grid points
# n            the number of sample to be generated in each data set
# weight       the proportion of samples to be generated for the mixture component corresponding to the minority of samples
# variances    the different variances to be simulated, specified as factor with respect to the MAD of the data
# GE_norm      which of the GE grid points should be assumed as the mean for the GE control data
# CN_norm      which of the CN inner grid points should be assumed as the mean for the CN control data
# seed         the seed
# call_probs   call probabilities to be specified for intCNGEan

simulation <- function(GE_data, CN_data, Inner=3:5, Outer=c(1:2,6:7), probs_GE=c(0.025,0.075,0.3,0.5,0.7,0.925,0.975), probs_CN=c(0.025,0.075,0.3,0.5,0.7,0.925,0.975), n=100, weight=1/10, variances=c(1/4,1/2,1,2,4), GE_norm=4, CN_norm=2, seed=888, call_probs=c(0.001,0.005,0.01,0.0125,0.04,0.925,0.99)){

if(length(Inner)+length(Outer) != length(probs_CN)){stop("The lengths of Inner and Outer jointly must be equivalent to the length of probs_CN.")}
#GE
GE_Werte <- log2(quantile(2^(as.numeric(unlist(GE_data))), probs=probs_GE))
GE_normal <- GE_Werte[GE_norm]
## CN
CN_Werte_aussen <- log2(quantile(2^(as.numeric(unlist(CN_data))), probs=probs_CN[Outer]))
CN_Werte_innen <- log2(quantile(2^(as.numeric(unlist(CN_data))), probs=probs_CN[Inner]))
CN_normal <- CN_Werte_innen[CN_norm]

log2Mad <- function(x){
    log2(mad(2^x))
}

log2Med <- function(x){
    log2(median(2^x))
}

CN_mad <- apply(CN_data,1,mad)
CN_mad_med <- log2Med(CN_mad)
GE_mad <- apply(GE_data,1,mad)
GE_mad_med <- log2Med(GE_mad)

Simulation <- function(n,CN_aussen,CN_innen,GE,nr_aussen,nr_innen,weight_wo,CN_mad_med,GE_mad_med){
    cn_nr_innen <- nr_innen %% length(Inner)
    if(cn_nr_innen == 0){cn_nr_innen <- length(Inner)}
    ge_nr_innen <- ceiling(nr_innen/length(Inner))
    
    cn_nr_aussen <- nr_aussen %% length(Outer)
    if(cn_nr_aussen == 0){cn_nr_aussen <- length(Outer)}
    ge_nr_aussen <- ceiling(nr_aussen/length(Outer))
    
    verteile <- round(n*(1-weight))
    if(weight_wo == 1){
        x <- c(rnorm(round(n*weight),CN_innen[cn_nr_innen],CN_mad_med),rnorm(round(n*(1-weight)),CN_aussen[cn_nr_aussen],CN_mad_med))
        y <- c(rnorm(round(n*weight),GE[ge_nr_innen],GE_mad_med),rnorm(round(n*(1-weight)),GE[ge_nr_aussen],GE_mad_med))
        }
    if(weight_wo == 2){  
        x <- c(rnorm(round(n*weight)/2,CN_innen[cn_nr_innen],CN_mad_med),rnorm(round(n*weight),CN_aussen[cn_nr_aussen],CN_mad_med),rnorm(round(n*(1-weight))-round(n*weight)/2,CN_innen[cn_nr_innen],CN_mad_med))
        y <- c(rnorm(round(n*weight)/2,GE[ge_nr_innen],GE_mad_med),rnorm(round(n*weight),GE[ge_nr_aussen],GE_mad_med),rnorm(round(n*(1-weight))-round(n*weight)/2,GE[ge_nr_innen],GE_mad_med))
    }
    x_r <- rnorm(2*n,CN_normal,CN_mad_med)  
    y_r <- rnorm(2*n,GE_normal,GE_mad_med) 
    
    output <- list(x=x,y=y,x_r=x_r,y_r=y_r)
    output
}

nr_aussen <- length(probs_GE)*length(Outer)
nr_innen <- length(probs_GE)*length(Inner)
nr_weights <- 2
nr_variance <- length(variances)

GE_ref_med <- GE_normal
CN_ref_med <- CN_normal

simulation_x <- array(numeric(n*nr_aussen*nr_innen*nr_weights*nr_variance), dim=c(n,nr_aussen,nr_innen,nr_weights,nr_variance))
simulation_y <- array(numeric(n*nr_aussen*nr_innen*nr_weights*nr_variance), dim=c(n,nr_aussen,nr_innen,nr_weights,nr_variance))

simulation_x_r <- array(numeric(2*n*nr_aussen*nr_innen*nr_weights*nr_variance), dim=c(2*n,nr_aussen,nr_innen,nr_weights,nr_variance))
simulation_y_r <- array(numeric(2*n*nr_aussen*nr_innen*nr_weights*nr_variance), dim=c(2*n,nr_aussen,nr_innen,nr_weights,nr_variance))
set.seed(seed)
z <- 1
for(a in 1:nr_aussen){
    for(i in 1:nr_innen){
        for(wo in 1:nr_weights){ 
            for(f in 1:nr_variance){
                    sim <- Simulation(n=n,CN_aussen=CN_Werte_aussen,CN_innen=CN_Werte_innen,GE=GE_Werte,nr_aussen=a,nr_innen=i,weight_wo=wo,CN_mad_med=variances[f]*CN_mad_med,GE_mad_med=variances[f]*GE_mad_med)#} #c(1/4,1/2,1,2,4)

                    simulation_x[,a,i,wo,f] <- sim$x
                    simulation_y[,a,i,wo,f] <- sim$y
                    
                    simulation_x_r[,a,i,wo,f] <- sim$x_r
                    simulation_y_r[,a,i,wo,f] <- sim$y_r
                    
                    z <- z + 1
             }
         }
     }
}

##################################################### DATENSATZ-GENERIERUNG ##############################################################
#add <- seq(length(probs_CN),length(probs_CN)*(length(probs_GE)-1),length(probs_CN))
#aussen <- c(Outer,Outer+add[1],Outer+add[2],Outer+add[3],Outer+add[4],Outer+add[5],Outer+add[6])
#innen <- c(Inner,Inner+add[1],Inner+add[2],Inner+add[3],Inner+add[4],Inner+add[5],Inner+add[6])

### Datensatz in "normaler" Datensatzform
simulation_x_geordnet <- matrix(numeric((nr_aussen*nr_innen*nr_weights*nr_variance)*n), ncol=n)
simulation_y_geordnet <- simulation_x_geordnet
simulation_x_r_geordnet <-  matrix(numeric((nr_aussen*nr_innen*nr_weights*nr_variance)*2*n), ncol=2*n)
simulation_y_r_geordnet <- simulation_x_r_geordnet


kombi <- matrix(numeric((nr_aussen*nr_innen*nr_weights*nr_variance)*4), ncol=4)
z <- 1
for(a in 1:nr_aussen){
     for(i in 1:nr_innen){
        for(wo in 1:nr_weights){ 
             for(f in 1:nr_variance){
                 kombi[z,] <- c(a,i,wo,f)
                 simulation_x_geordnet[z,] <- simulation_x[,a,i,wo,f]
                 simulation_y_geordnet[z,] <- simulation_y[,a,i,wo,f]
                 simulation_x_r_geordnet[z,] <- simulation_x_r[,a,i,wo,f]
                 simulation_y_r_geordnet[z,] <- simulation_y_r[,a,i,wo,f]
                 z <- z + 1
             }
         }    
     }
}

rownames(simulation_x_geordnet) <- rownames(simulation_y_geordnet) <- rownames(simulation_x_r_geordnet) <- rownames(simulation_y_r_geordnet) <- as.character(1:nrow(simulation_x_geordnet))
#z1 <- z-1
################################################################ x #########################################################################
# Einmal
simulation_x_doppelt_geordnet1 <- matrix(numeric(nrow(simulation_x_geordnet)*ncol(simulation_x_geordnet)), ncol=ncol(simulation_x_geordnet))
simulation_y_doppelt_geordnet1 <- matrix(numeric(nrow(simulation_y_geordnet)*ncol(simulation_y_geordnet)), ncol=ncol(simulation_y_geordnet))
simulation_x_r_doppelt_geordnet1 <- matrix(numeric(nrow(simulation_x_r_geordnet)*ncol(simulation_x_r_geordnet)), ncol=ncol(simulation_x_r_geordnet))
simulation_y_r_doppelt_geordnet1 <- matrix(numeric(nrow(simulation_x_r_geordnet)*ncol(simulation_x_r_geordnet)), ncol=ncol(simulation_x_r_geordnet))

sortiert <- sort(apply(simulation_x_geordnet,1,median), index.return=TRUE)

for(j in 1:ncol(simulation_x_geordnet)){
#print(j)
     simulation_x_doppelt_geordnet1[,j] <- simulation_x_geordnet[sortiert$ix,j]
     simulation_y_doppelt_geordnet1[,j] <- simulation_y_geordnet[sortiert$ix,j]
     #simulation_x_r_doppelt_geordnet1[,j] <- simulation_x_r_geordnet[,j]
     #simulation_y_r_doppelt_geordnet1[,j] <- simulation_y_r_geordnet[,j]
}

for(j in 1:ncol(simulation_x_r_geordnet)){
     simulation_x_r_doppelt_geordnet1[,j] <- simulation_x_r_geordnet[,j]
     simulation_y_r_doppelt_geordnet1[,j] <- simulation_y_r_geordnet[,j]
}

rownames(simulation_x_doppelt_geordnet1 ) <- rownames(simulation_y_doppelt_geordnet1) <- rownames(simulation_x_r_doppelt_geordnet1) <- rownames(simulation_y_r_doppelt_geordnet1) <- rownames(simulation_x_geordnet)[sortiert$ix]

ind <- 1:nrow(simulation_x_doppelt_geordnet1)
teile <- 16
breaks <- c(0,(1:teile)*length(ind)/teile)
for(i in 1:(length(breaks)-1)){
   #ind[(breaks[i]+1):breaks[i+1]] <- sample(ind[(breaks[i]+1):breaks[i+1]],length(ind[(breaks[i]+1):breaks[i+1]]))
   sortiert_breaks <- sort(apply(simulation_y_doppelt_geordnet1[(breaks[i]+1):breaks[i+1],],1,median), index.return=TRUE)
   ind[(breaks[i]+1):breaks[i+1]] <- ind[(breaks[i]+1):breaks[i+1]][sortiert_breaks$ix]
}
ind_rueck <- sort(ind, index.return=TRUE)

#Plot
#k <- 6
#zufall_patient_x <- simulation_x_doppelt_geordnet1[ind,k]
#plot(1:nrow(simulation_x_geordnet),zufall_patient_x)

# Daten ordnen
simulation_x_doppelt_geordnet <- simulation_x_doppelt_geordnet1[ind,]
simulation_y_doppelt_geordnet <- simulation_y_doppelt_geordnet1[ind,]
simulation_x_r_doppelt_geordnet <- simulation_x_r_doppelt_geordnet1
simulation_y_r_doppelt_geordnet <- simulation_y_r_doppelt_geordnet1
rownames(simulation_x_doppelt_geordnet) <- rownames(simulation_y_doppelt_geordnet) <- rownames(simulation_x_r_doppelt_geordnet) <- rownames(simulation_y_r_doppelt_geordnet) <- rownames(simulation_x_doppelt_geordnet1)[ind]

#########################################################################################################################
################################################# Ratios bilden #########################################################
cn_norm_median <- apply(simulation_x_r_doppelt_geordnet[,1:n],1,median)
simulation_x_doppelt_geordnet_ratios <- simulation_x_doppelt_geordnet - cn_norm_median
simulation_x_r_doppelt_geordnet_ratios <- simulation_x_r_doppelt_geordnet[,(n+1):(2*n)] - cn_norm_median

ge_norm_median <- apply(simulation_y_r_doppelt_geordnet[,1:n],1,median)
simulation_y_doppelt_geordnet_ratios <- simulation_y_doppelt_geordnet - ge_norm_median
simulation_y_r_doppelt_geordnet_ratios <- simulation_y_r_doppelt_geordnet[,(n+1):(2*n)] - ge_norm_median

#########################################################################################################################
#########################################################################################################################
####################################################### Segmentierung mit CBS ###########################################
#########################################################################################################################
#########################################################################################################################
require(DNAcopy)
CNA <- CNA(genomdat=simulation_x_doppelt_geordnet_ratios, chrom=rep(1,nrow(simulation_x_doppelt_geordnet_ratios)), maploc=1:nrow(simulation_x_doppelt_geordnet_ratios), data.type="logratio",sampleid=1:n)
smoothed.CNA <- smooth.CNA(CNA)
segment.smoothed.CNA <- DNAcopy::segment(smoothed.CNA, verbose = 1)

samples_liste <- strsplit(segment.smoothed.CNA$output$ID, split="X")
nehme_2 <- function(x){
    x[2]
}
samples <- as.numeric(lapply(samples_liste, nehme_2))

simulation_x_doppelt_geordnet_ratios_segmentiert <- simulation_x_doppelt_geordnet_ratios
for(s in 1:ncol(simulation_x_doppelt_geordnet_ratios_segmentiert)){
    welche <- which(samples == s)
    for(w in welche){
        spanne <- segment.smoothed.CNA$output$loc.start[w]:segment.smoothed.CNA$output$loc.end[w]
        simulation_x_doppelt_geordnet_ratios_segmentiert[spanne,s] <- segment.smoothed.CNA$output$seg.mean[w] #median(simulation_x_doppelt_geordnet_segmentiert[spanne,s])
    }
}

###############################

CNA_r <- CNA(genomdat=simulation_x_r_doppelt_geordnet_ratios, chrom=rep(1,nrow(simulation_x_r_doppelt_geordnet_ratios)), maploc=1:nrow(simulation_x_r_doppelt_geordnet_ratios), data.type="logratio",sampleid=1:n)
smoothed.CNA_r <- smooth.CNA(CNA_r)
segment.smoothed.CNA_r <- DNAcopy::segment(smoothed.CNA_r, verbose = 1)

samples_r_liste <- strsplit(segment.smoothed.CNA_r$output$ID, split="X")
nehme_2 <- function(x){
    x[2]
}
samples_r <- as.numeric(lapply(samples_r_liste, nehme_2))

simulation_x_r_doppelt_geordnet_ratios_segmentiert <- simulation_x_r_doppelt_geordnet_ratios
for(s in 1:ncol(simulation_x_r_doppelt_geordnet_ratios_segmentiert)){
    welche <- which(samples_r == s)
    for(w in welche){
        spanne_r <- segment.smoothed.CNA_r$output$loc.start[w]:segment.smoothed.CNA_r$output$loc.end[w]
        simulation_x_r_doppelt_geordnet_ratios_segmentiert[spanne_r,s] <- segment.smoothed.CNA_r$output$seg.mean[w] #<- median(simulation_x_r_doppelt_geordnet_segmentiert[spanne_r,s])
    }
}

####################################### artificial call probabilities, e.g., fot intCNGEan ########################################
rueck1 <- sort(sortiert$ix, index.return=TRUE)
probs <- call_probs
probs_norm <- 1-call_probs-call_probs[length(call_probs):1]

zahlen_umgeordnet <- 1:nrow(kombi)
zahlen_normal <- zahlen_umgeordnet[ind_rueck$ix][rueck1$ix]
sortierung <- sort(zahlen_normal, index.return=TRUE)

gainprob <- function(zeile){
    a <- kombi[sortierung$ix,][zeile,1]
    i <- kombi[sortierung$ix,][zeile,2]
    wo <- kombi[sortierung$ix,][zeile,3]
    
    cn_nr_innen <- i %% 3
    if(cn_nr_innen == 0){cn_nr_innen <- 3}
    
    cn_nr_aussen <- a %% 4
    if(cn_nr_aussen == 0){cn_nr_aussen <- 4}
        
    aussen <- probs[c(1,2,6,7)]
    innen <- probs[3:5]

    if(wo == 1){out <- c(rep(innen[cn_nr_innen],round(n*weight)), rep(aussen[cn_nr_aussen],n*(1-weight))) }
    if(wo == 2){out <- c(rep(innen[cn_nr_innen],round(n*weight)/2), rep(aussen[cn_nr_aussen],round(n*weight)), rep(innen[cn_nr_innen],round(n*(1-weight))-round(n*weight)/2)) }
    out
}

lossprob <- function(zeile){
    a <- kombi[sortierung$ix,][zeile,1]
    i <- kombi[sortierung$ix,][zeile,2]
    wo <- kombi[sortierung$ix,][zeile,3]
    
    cn_nr_innen <- i %% 3
    if(cn_nr_innen == 0){cn_nr_innen <- 3}
    
    cn_nr_aussen <- a %% 4
    if(cn_nr_aussen == 0){cn_nr_aussen <- 4}
       
    aussen <- probs[c(7,6,2,1)]
    innen <- probs[5:3]
    
    if(wo == 1){out <- c(rep(innen[cn_nr_innen],round(n*weight)), rep(aussen[cn_nr_aussen],n*(1-weight))) }
    if(wo == 2){out <- c(rep(innen[cn_nr_innen],round(n*weight)/2), rep(aussen[cn_nr_aussen],round(n*weight)), rep(innen[cn_nr_innen],round(n*(1-weight))-round(n*weight)/2)) }
    out
}

normprob <- function(zeile){
    a <- kombi[sortierung$ix,][zeile,1]
    i <- kombi[sortierung$ix,][zeile,2]
    wo <- kombi[sortierung$ix,][zeile,3]
    
    cn_nr_innen <- i %% 3
    if(cn_nr_innen == 0){cn_nr_innen <- 3}
    
    cn_nr_aussen <- a %% 4
    if(cn_nr_aussen == 0){cn_nr_aussen <- 4}
        
    aussen <- probs_norm[c(1,2,6,7)]
    innen <- probs_norm[3:5]
    
    if(wo == 1){out <- c(rep(innen[cn_nr_innen],round(n*weight)), rep(aussen[cn_nr_aussen],n*(1-weight))) }
    if(wo == 2){out <- c(rep(innen[cn_nr_innen],round(n*weight)/2), rep(aussen[cn_nr_aussen],round(n*weight)), rep(innen[cn_nr_innen],round(n*(1-weight))-round(n*weight)/2)) }
    out
}

gainprobs <- matrix(numeric(nrow(kombi)*n), nrow=nrow(kombi))
lossprobs <- matrix(numeric(nrow(kombi)*n), nrow=nrow(kombi))
normprobs <- matrix(numeric(nrow(kombi)*n), nrow=nrow(kombi))

for(j in 1:nrow(kombi)){
    gainprobs[j,] <- gainprob(j)
}

for(j in 1:nrow(kombi)){
    lossprobs[j,] <- lossprob(j)
}

for(j in 1:nrow(kombi)){
    normprobs[j,] <- normprob(j)
}

call_probs <- list(gain=gainprobs, loss=lossprobs,norm=normprobs)

############################################### Output erstellen ##################################################################
results <- list(CN=simulation_x_doppelt_geordnet_ratios,CN.norm=simulation_x_r_doppelt_geordnet_ratios,CN.seg=simulation_x_doppelt_geordnet_ratios_segmentiert,CN.norm.seg=simulation_x_r_doppelt_geordnet_ratios_segmentiert,GE=simulation_y_doppelt_geordnet_ratios,GE.norm=simulation_y_r_doppelt_geordnet_ratios,callprobs=call_probs)
results
}
