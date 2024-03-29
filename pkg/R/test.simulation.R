# This script is part of the intcomp project: 
# http://intcomp.r-forge.r-project.org/
# License: FreeBSD, http://en.wikipedia.org/wiki/BSD_licenses
# Copyright 2011 Leo Lahti and Martin Schafer, <leo.lahti@iki.fi>. All
# rights reserved.


# GE           either a matrix of gene expression data, without annotation, or a list containing an object 'data' (used only for method="schaefer"). If nothing is supplied, data from Schaefer et al. (2009) is taken.
# CN           either a matrix of copy number data, without annotation, or a list containing an object 'data' (used only for method="schaefer"). If nothing is supplied, data from Schaefer et al. (2009) is taken.
# probespanGE  the width of gene expression probes
# probespanCN  the width of copy number probes
# method       either "schaefer" or "ferrari". The option "schaefer" performs simulation following the approach presented in Schaefer et al. (2009): Integrated analysis of copy number alterations and gene expression: a bivariate assessment of equally directed abnormalities. Bioinformatics 25(24):3228�3235.
#              The option "ferrari" leads to use of data simulated by Francesco Ferrari.
# Inner        indices of inner grid points for simulation of the copy number data (used only for method="schaefer")
# Outer        indices of outer grid points for simulation of the copy number data (used only for method="schaefer")
#              (one mixture component of the normal mixture distribution has as copy number coordinate of its mean an inner grid point, the other mixture component has as copy number coordinate of its mean an outer grid point)
# probs_GE     the quantiles used to calculate the GE grid points (used only for method="schaefer")
# probs_CN     the quantiles used to calculate the CN grid points (used only for method="schaefer")
# cancer_GE    specifies the gene expression coordinates of the points in the quantile grid which should be defined as cancer genes (used only for method="schaefer")
# cancer_CN    specifies the copy number coordinates of the points in the quantile grid which should be defined as cancer genes (used only for method="schaefer")
# n            the number of sample to be generated in each data set (used only for method="schaefer")
# weight       the proportion of samples to be generated for the mixture component corresponding to the minority of samples (used only for method="schaefer")
# variances    the different variances to be simulated, specified as factor with respect to the MAD of the data (used only for method="schaefer")
# GE_norm      which of the GE grid points should be assumed as the mean for the GE control data (used only for method="schaefer")
# CN_norm      which of the CN inner grid points should be assumed as the mean for the CN control data (used only for method="schaefer")
# seed         the seed (used only for method="schaefer")
# call_probs   call probabilities to be specified for intCNGEan (used only for method="schaefer")

test.simulation <- function(GE=NULL, CN=NULL, probespanGE = 16, probespanCN = 16, method, Inner=3:5, Outer=c(1:2,6:7),
probs_GE=c(0.025,0.075,0.3,0.5,0.7,0.925,0.975), probs_CN=c(0.025,0.075,0.3,0.5,0.7,0.925,0.975),
cancer_GE=c(1,1,2,2,6,6,7,7),cancer_CN=c(1,2,1,2,3,4,3,4), n=100, weight=1/10, variances=c(1/4,1/2,1,2,4),
GE_norm=4, CN_norm=2, seed=42, call_probs=c(0.001,0.005,0.01,0.0125,0.04,0.925,0.99)){

if(method=="schaefer"){
if(length(GE) == 0 | length(CN) == 0){
require(ediraAMLdata)
data(AMLdata, package="ediraAMLdata")
}
if(is.list(GE)){
    if(length(GE$data) > 0){GE <- as.matrix(GE$data)} else {stop("If GE is a list, it must contain an element 'data'.")}
}
if(is.list(CN)){
    if(length(CN$data) > 0){CN <- as.matrix(CN$data)} else {stop("If CN is a list, it must contain an element 'data'.")}
}
if(is.matrix(GE) & is.matrix(CN)){
if(length(probs_GE) != length(probs_CN)){stop("probs_GE and probs_CN must have equal length.")}
if(length(cancer_GE) != length(cancer_CN)){stop("cancer_GE and cancer_CN must have equal length.")}

sim <- simulation(GE, CN, Inner, Outer, probs_GE, probs_CN, n, weight, variances, GE_norm, CN_norm, seed, call_probs)
################################################################################################################

# annotations for simulated data
mittelpunkte <- (1:nrow(sim$GE))*100
ge_anfang <- mittelpunkte-probespanGE
cn_anfang <- mittelpunkte-probespanCN
ge_ende <- mittelpunkte+probespanGE
cn_ende <- mittelpunkte+probespanGE
chr <- rep(1,nrow(sim$GE))
ge_info <- data.frame(start=ge_anfang,end=ge_ende, chr=chr, loc=mittelpunkte)
cn_info <- data.frame(start=cn_anfang,end=cn_ende, chr=chr, loc=mittelpunkte)
rownames(ge_info) <- rownames(sim$GE)

# patient data with copy number data segmented
ge <- list(data=sim$GE, info=ge_info) 
cn <- list(data=sim$CN.seg, info=cn_info)

# patient and reference data with copy number data segmented
ge.norm <- list(data=sim$GE.norm, info=ge_info)
cn.norm <- list(data=sim$CN.norm.seg, info=cn_info)  

# raw (unsegmented copy number data)
cn.raw <- list(data=sim$CN, info=cn_info)

# add column names
colnames(ge$data) <- as.character(1:ncol(ge$data))
colnames(cn$data) <- as.character(1:ncol(cn$data))
colnames(cn.raw$data) <- as.character(1:ncol(cn.raw$data))

# Labels
Labels <- c(rep(1,100),rep(0,100))

# define cancer genes
matrix_dimension <- length(probs_GE)
nr_scenarios <- 2*length(variances)
nr_outer <- length(Outer)*matrix_dimension
nr_inner <- length(Inner)*matrix_dimension
anfang <- (seq(0,nr_inner*nr_scenarios*nr_outer,nr_inner*nr_scenarios)+1)
ende <- seq(0,nr_inner*nr_scenarios*nr_outer,nr_inner*nr_scenarios)[-1]
a <- c(1,1+seq(nr_scenarios,nr_inner*nr_scenarios*length(cancer_GE),nr_scenarios))
addiere <- function(x){a + x}
x <- 0:(length(variances)-1)
alles <- sort(as.numeric(sapply(x, addiere)))
#alles <- sort(a+seq(0,length(variances)-1,1)) #c(a,a+1,a+2,a+3,a+4)

gib_gesamt_index <- function(index_ge, index_cn){(index_ge-1)*length(Outer)+index_cn}

cancerGenes_temp <- numeric(0)
for(i in 1:length(cancer_GE)){
        index <- gib_gesamt_index(cancer_GE[i],cancer_CN[i])
        cancerGenes_temp <- c(cancerGenes_temp, anfang[index]:ende[index])   
}
cancerGenes <- as.character(cancerGenes_temp)[alles][1:(nr_inner*length(variances)*length(cancer_GE))]

# call CN data
CN <- process.copynumber(cn.raw, cn.seg = NULL, probespanCN = 16, prior = "all", organism = "human")
cn.call <- list(data=assayDataElement(CN, 'calls'), info=cn_info)

out=list(ge=ge, cn=cn, ge.norm=ge.norm, cn.norm=cn.norm, cn.raw=cn.raw, cn.call=cn.call, cn.cghCall=CN, Labels=Labels, cancerGenes=cancerGenes, callprobs=sim$callprobs)
}
}
####################################################################################### Ferrari simulations ####################################################################################
if(method=="ferrari"){
data(ferrari_simulations)
# annotations for simulated data
ge_data <- as.matrix(GE$data)
ge_norm_data <- as.matrix(GE.norm$data)
cn_data <- as.matrix(CN.seg$data)
cn_norm_data <- as.matrix(CN.norm.seg$data)
cn_raw_data <- as.matrix(CN$data)
rownames(ge_data) <- as.character(1:nrow(ge_data))
rownames(ge_norm_data) <- as.character(1:nrow(ge_norm_data))    
rownames(cn_data) <- as.character(1:nrow(cn_data))
rownames(cn_norm_data) <- as.character(1:nrow(cn_norm_data))
rownames(cn_raw_data) <- as.character(1:nrow(cn_raw_data))

mittelpunkte <- GE$info$position
anfang <- mittelpunkte-probespanGE
ende <- mittelpunkte+probespanGE
chr <- GE$info$chromosome
ge_info <- data.frame(start=anfang,end=ende, chr=chr, loc=mittelpunkte)

mittelpunkte <- CN.seg$info$position
anfang <- mittelpunkte-probespanCN
ende <- mittelpunkte+probespanCN
chr <- CN.seg$info$chromosome
cn_info <- data.frame(start=anfang,end=ende, chr=chr, loc=mittelpunkte)

mittelpunkte <- CN$info$position
anfang <- mittelpunkte-probespanCN
ende <- mittelpunkte+probespanCN
chr <- CN$info$chromosome
cn_raw_info <- data.frame(start=anfang,end=ende, chr=chr, loc=mittelpunkte)
rownames(ge_info) <- rownames(ge_data)
rownames(cn_info) <- rownames(cn_data)
rownames(cn_raw_info) <- rownames(cn_raw_data)

# patient data with copy number data segmented
ge <- list(data=ge_data, info=ge_info)
cn <- list(data=cn_data, info=cn_info)

# patient and reference data with copy number data unsegmented
ge.norm <- list(data=ge_norm_data, info=ge_info)
cn.norm <- list(data=cn_norm_data, info=cn_info)

# raw (unsegmented copy number data)
cn.raw <- list(data=cn_raw_data, info=cn_raw_info)

# Labels
Labels <- c(rep(1,6),rep(0,6))

################################
# define cancer genes
chr <- c(1,1,3,3,5,6,10,17,18,19)
starts <- c(20,150,30,120,135,30,20,40,60,40)

List_all <- list()
List <- list()
for(j in 1:length(starts)){
    List[[j]] <- which(GE$info$chromosome == chr[j] & GE$info$position >= starts[j]*10^6 & GE$info$position <= (starts[j]+10)*10^6)
    List_all <- c(List_all, List[[j]])
}

cancerGenes <- as.character(List_all)

# call CN data
CN <- process.copynumber(cn.raw, cn.seg = NULL, 
                         probespanCN = 100, 
			 prior = "all", organism = "human")

cn.call <- list(data=assayDataElement(CN, 'calls'), info = cn.raw$info)

# For each ge probe, find the closest cn.call segment
cn.inds <- c()
for (i in 1:nrow(ge$info)) {
  chr <- ge$info[i, "chr"]
  loc <- ge$info[i, "start"]
  inds <- which(cn.call$info$chr == chr) # For each ge probe, find the closest cn.call segment
  starts <- cn.call$info[inds, "start"]
  ind <- which.min(abs(starts - loc))
  cn.inds[[i]] <- ind
}
cn.call.matched <- list()
cn.call.matched$data <- as.matrix(cn.call$data)[cn.inds,]
cn.call.matched$info <- ge$info

out <- list(ge = ge, 
	    cn.raw = cn.raw,
            cn.seg = cn, 	    
	    cn.call = cn.call.matched,
	    cn.cghCall = CN, 
	    ge.norm = ge.norm, 
	    cn.norm = cn.norm, 
	    Labels = Labels, 
	    cancerGenes = cancerGenes)

}

  return(out)

}
