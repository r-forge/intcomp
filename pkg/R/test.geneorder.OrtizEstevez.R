test.geneorder.OrtizEstevez <- function (ge, cn) {

  require("DNAcopy")

  ordg <- getOrderedGenes(cn, ge, N = nrow(ge$data))

  ordg
}


getOrderedGenes <- function(cn, ge, N=300)
{	
	umbral_pos <- 0.10
	umbral_neg <- -0.10
	
	allDataGE <- list();
	for(ii in 1:(dim(ge$data)[2])){
 		aux <- CNA(ge$data[,ii], as.numeric(ge$info$chr), as.numeric(ge$info$loc), data.type="logratio");
 		auxSeg <- segment(aux)
 		allDataGE <- c(allDataGE,list(auxSeg));
	}

	sGE <- matrixFromRegions(allDataGE);
	sGEn <- apply(sGE, 2, function(x) x-median(x))

	allDataCN <- list();
	for(ii in 1:(dim(cn$data)[2])){
 		aux <- CNA(cn$data[,ii], as.numeric(cn$info$chr), as.numeric(cn$info$loc), data.type="logratio");
 		auxSeg <- segment(aux)
		allDataCN <- c(allDataCN,list(auxSeg));
	}

	sCN <- matrixFromRegions(allDataCN);
	sCNn <- apply(sCN, 2, function(x) x-median(x))

	stateCN <- apply(sCNn, 2, function(x){ 
							x[which(x >= umbral_pos)] <- 1
							x[which(x <= umbral_neg)] <- -1
							x[which(x > umbral_neg & x < umbral_pos)] <- 0
							return(x)
						})
								
	stateGE <- apply(sGEn, 2, function(x){ 
							x[which(x >= umbral_pos)] <- 1
							x[which(x <= umbral_neg)] <- -1
							x[which(x > umbral_neg & x < umbral_pos)] <- 0
							return(x)
						})


	gainCN <- matrix(0,nrow(sCNn),ncol(sCNn))
	rownames(gainCN) <- rownames(cn$info)
	gainCN[stateCN == 1] <- sCNn[stateCN == 1]
	gainCN <- apply(gainCN, 1, sum)

	lossCN <- matrix(0,nrow(sCNn),ncol(sCNn))
	rownames(lossCN) <- rownames(cn$info)
	lossCN[stateCN == -1] <- sCNn[stateCN == -1]
	lossCN <- apply(lossCN, 1, sum)

	gainGE <- matrix(0,nrow(sGEn),ncol(sGEn))
	rownames(gainGE) <- rownames(ge$info)
	gainGE[stateGE == 1] <- sGEn[stateGE == 1]
	gainGE <- apply(gainGE, 1, sum)

	lossGE <- matrix(0,nrow(sGEn),ncol(sGEn))
	rownames(lossGE) <- rownames(ge$info)
	lossGE[stateGE == -1] <- sGEn[stateGE == -1]
	lossGE <- apply(lossGE, 1, sum)

	gainlossCNgenes <- getRanking(apply(cbind(gainCN,abs(lossCN)), 1, max))
	gainlossGEgenes <- getRanking(apply(cbind(gainGE,abs(lossGE)), 1, max))

	orderCN <- 1:length(gainlossCNgenes)
	names(orderCN) <- gainlossCNgenes

	orderGE <- 1:length(gainlossGEgenes)
	names(orderGE) <- gainlossGEgenes

	allgenes <- cbind(CN=orderCN[names(orderCN)], GE=orderGE[names(orderCN)])
	rownames(allgenes) <- names(orderCN)
	allgenes <- cbind(allgenes, sum=apply(allgenes, 1, sum))

	allgenesOrdered <- names(sort(allgenes[,3]))[1:N]
	
	return (allgenesOrdered)
}


matrixFromRegions <- function(cbsData){
  nSamples <- length(cbsData);
  
  nbrLoci <- length(cbsData[[1]]$data$chrom);
  dataMat <- matrix(data=0, nrow = nbrLoci, ncol = nSamples);

  nChr <- sort(unique(as.numeric(cbsData[[1]]$data$chrom)));
  for (i in 1:nSamples){
    sData <- cbsData[[i]];
    lociDone <- 0;
    for (j in 1:length(nChr)){
      indChrom <- sData$data$chrom==j;
      lociPos <- sData$data$maploc[indChrom];
      
      indChrLoci <- sData$output$chrom==j;
      scData <- sData$output[indChrLoci,];

      indReg <- 1;
      nLoci <- length(lociPos);
      auxMat <- matrix(data=0,nrow=nLoci,ncol=1);
      for (k in 1:nLoci){
        locusPos <- lociPos[k];
        found = FALSE;
        while(indReg < sum(indChrLoci) & found == FALSE){
          if(locusPos < scData$loc.start[indReg] & (indReg ==1 || (indReg>1 & locusPos > scData$loc.end[indReg-1]))){
            auxMat[k,1] <- scData$seg.mean[indReg];
            found = TRUE;    
          }else{
            if(locusPos >= scData$loc.start[indReg] & locusPos <= scData$loc.end[indReg]){
              auxMat[k,1] <- scData$seg.mean[indReg];
              found = TRUE;
            }else{
              indReg <- indReg + 1;
            }
          }
        }
        if (found == FALSE & indReg == sum(indChrLoci)){
          auxMat[k,1] <- scData$seg.mean[indReg];
        }
      }
      dataMat[(lociDone+1):(lociDone+nLoci),i]  <- auxMat;                   
      lociDone <- lociDone+nLoci;
    }
  }
  return(dataMat);
}

getRanking <- function(x){
	rankedGenes <- c()
	while(max(x) != 0){
		rankedGenes <- c(rankedGenes, names(which(x == max(x))))
		x[which(x == max(x))] <- 0 
	}
	return (rankedGenes)
}





