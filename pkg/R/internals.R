
removeAFFX <- function (sets) {                                      
  # Remove AFFX control sets from set list            
  sets[!findAFFX(sets)]                
}
                                                              
                           
findAFFX <- function (sets) {              
  # Find AFFX control sets from the given set list                                                        
  affx.set<-logical(length(sets))    
  for (i in 1:length(sets)) {                
    affx.set[[i]]<-(substr(sets[[i]],1,4)=="AFFX")                      
  }                                      
  affx.set                                              	        
}    



get.pma.table <- function (cn, ge, typex, typez) {

   # Define data matrices
   dna <- t(cn$data)
   rna <- t(ge$data)
	             
   # Find optimal model parameters
   perm.out <- CCA.permute(x = rna,
			   z = dna,
                           typex = typex,
                           typez = typez)
			                          
    # Calculate model
    out <- CCA(x = rna,
               z = dna,
               typex = typex,
	       typez = typez,
	    penaltyx = perm.out$bestpenaltyx,
		   v = perm.out$v.init,
	    penaltyz = perm.out$bestpenaltyz)
                     
     tab <- as.data.frame(list(rownames(cn$data),
			  cn$info$chr,
			  cn$info$loc,
			 score = I(as.numeric(abs(out$u)))))
             # Selecting this score (out$u) based on the example in ?CCA
             # taking absolute score since both gain/loss can be observed
             # and ground truth list does not distinguish between these.
     colnames(tab) <- c("probeid", "chr", "loc", "score")
     rownames(tab) <- rownames(cn$data)

     tab
	         
}
		     
		     


centerData <- function (X, rm.na = FALSE, meanvalue = NULL) {

  # Center columns of a data matrix to desired value(s).
  # Estimate parameters (mean and variance) separately for each column.
  # (C) Leo Lahti 2010. License: FreeBSD.

  if (!rm.na) {
    xcenter <- colMeans(X)
    X2 <- X - rep(xcenter, rep.int(nrow(X), ncol(X)))
  } else {  
    X2 <- array(NA, dim = c(nrow(X), ncol(X)), dimnames = dimnames(X))
    for (i in 1:ncol(X)) {
      x <- X[,i]
      nainds <- is.na(x)
      xmean <- mean(x[!nainds])
      X2[!nainds,i] <- x[!nainds] - xmean   
    }
    dimnames(X2) <- dimnames(X)
  }

  if (!is.null(meanvalue)) {
    # Shift the data to a specified value
    X2 <- X2 + meanvalue
  }


  
  X2
}

closest <- function(a, vec){which.min(abs(a - vec))}

esort <- function(x, sortvar, ...) {
  #Sort data frame dd by columns like: esort(dd, -z, b)

  attach(x)
  x <- x[with(x,order(sortvar,...)),]
  return(x)
  detach(x)
}

form2 <- function(X, samplenames, rm.chr = NULL){
  chr <- as.character(X$info$chr)
  chr[chr == 'X'] <- 23
  chr[chr == 'Y'] <- 24
  chr <- as.numeric(chr)

  remove.probes <- (chr %in% rm.chr)
  if (sum(remove.probes) > 0) {
    X$info <- X$info[!remove.probes,]
    X$data <- X$data[!remove.probes,]
    chr <- chr[!remove.probes]
  }

  ret <- cbind(as.data.frame(X$data[,colnames(X$data) %in% samplenames]), chr, X$info$loc)
  colnames(ret) <- c(colnames(ret)[1:(length(colnames(ret))-2)], "chromosome", "position")
  return(ret)
}
  

get.entrez.info <- function ( geneids ) {

  # Get Entrez GeneID chromosomal location information from entrez
  require(biomaRt)
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  geneid.info <- getBM(attributes = c("entrezgene", "start_position", "end_position", "chromosome_name"), values = geneids, mart = ensembl)

  # Remove duplicated geneids
  genes.uniq <- names(which(table(geneid.info[, "entrezgene"]) == 1))
  #geneid.info <- subset(geneid.info, entrezgene %in% genes.uniq)
  geneid.info <- geneid.info[geneid.info[["entrezgene"]] %in% genes.uniq, ]

  # Pick and arrange information
  info <- geneid.info[, c("start_position", "end_position", "chromosome_name")]
  colnames(info) <- c("start","end","chr")
  rownames(info) <- geneid.info[, "entrezgene"]

  info
  
}






get.neighboring.probes <- function (X, Y, chr, max.dist, control.arms = TRUE) {

  xinds <- yinds <- c()
  
  # Use arm information if it is available and not blocked
  if (("arm" %in% names(X$info)) && ("arm" %in% names(Y$info)) && control.arms) {    

    for (arm in c('p', 'q')){      
      # Investigate specified arm
      xchrinds <- which(as.character(X$info$chr) == chr & X$info$arm == arm)
      ychrinds <- which(as.character(Y$info$chr) == chr & Y$info$arm == arm)
      inds <- get.neighs(X, Y, xchrinds, ychrinds, max.dist)
      xinds <- c(xinds, inds$xinds)
      yinds <- c(yinds, inds$yinds)
    }
  } else {
    # Investigate the whole chromosome
    xchrinds <- which(as.character(X$info$chr) == chr)
    ychrinds <- which(as.character(Y$info$chr) == chr)
    inds <- get.neighs(X, Y, xchrinds, ychrinds, max.dist)
    xinds <- c(xinds, inds$xinds)
    yinds <- c(yinds, inds$yinds)
  }

  # return the indices
  list(xinds = xinds, yinds = yinds)
  
}


get.neighs <- function (X, Y, xchrinds, ychrinds, max.dist) {

  xinds <- yinds <- NULL

  if ( length(xchrinds) > 0 && length(ychrinds) > 0 ){
    
    #Find indices of closest probe from Y for each from X
    xi <- 1:length(xchrinds)
    yi <- sapply(as.numeric(as.character(X$info$loc[xchrinds])), closest, vec = as.numeric(as.character(Y$info$loc[ychrinds])))

    # Remove duplicates
    keep <- !duplicated(yi)
    xi <- xi[keep]
    yi <- yi[keep]
      
    # Corresponding indices between X and Y
    xinds <- xchrinds[xi]
    yinds <- ychrinds[yi]

    # delete indices which are further from each other than threshold
    near <- (abs(X$info$loc[xinds] - Y$info$loc[yinds]) < max.dist)
    xinds <- xinds[near]
    yinds <- yinds[near]
  
    # calculate mean location for each pair for ordering of the pairs
    xy.loc <- (X$info$loc[xinds] + Y$info$loc[yinds])/2

    # ensure probes are ordered by location
    ord <- order(xy.loc)
    xinds <- xinds[ord]
    yinds <- yinds[ord]

  }

  list(xinds = xinds, yinds = yinds)
  
}


impute <- function (X) {

  # Impute missing values from a Gaussian.
  # Estimate parameters (mean and variance) separately for each column.
  # (C) Leo Lahti 2010. License: FreeBSD.
  
  for (i in 1:ncol(X)) {
    x <- X[, i]
    nas <- is.na(x)
    X[nas, i] <- rnorm(sum(nas), mean(x[!is.na(x)]), sd(x[!is.na(x)]))
  }
  
  X
  
}

pairmatch <- function(X, Y, max.dist = 1e7, chrs = NULL){

  # Match probes between two data sets based on location information
  
  if (all(is.na(X$data))) {stop("X data is empty/NA.")}
  if (all(is.na(Y$data))) {stop("Y data is empty/NA.")}
  
  # If same number of rows and columns, assume that they match between data and info fields
  # if the names are completely non-overlapping (as in our example data set)
  if (nrow(X$data) == nrow(X$info) && length(intersect(rownames(X$info), rownames(X$data)))==0) {
    rownames(X$info) <- rownames(X$data)
  }
  if (nrow(Y$data) == nrow(Y$info) && length(intersect(rownames(Y$info), rownames(Y$data)))==0) {
    rownames(Y$info) <- rownames(Y$data)
  }

  # provide information in a form that has corresponding rows in
  # data and info matrices (only take those available in both)
  coms <- intersect(rownames(X$data), rownames(X$info))
  coms <- setdiff(coms, c(""))
  X$data <- X$data[coms,]
  X$info <- X$info[coms,]

  coms <- intersect(rownames(Y$data), rownames(Y$info))
  coms <- setdiff(coms, c(""))
  Y$data <- Y$data[coms,]
  Y$info <- Y$info[coms,]

  X$info[["chr"]] <- as.character(X$info[["chr"]])
  Y$info[["chr"]] <- as.character(Y$info[["chr"]])

  # X/Y chromosome for data set X
  if ("X" %in% X$info[["chr"]]) {X$info[X$info[["chr"]] == "X", "chr"] <- "23"}
  if ("Y" %in% X$info[["chr"]]) {X$info[X$info[["chr"]] == "Y", "chr"] <- "24"}

  # X/Y chromosome for data set X
  if ("X" %in% Y$info[["chr"]]) {Y$info[Y$info[["chr"]] == "X", "chr"] <- "23"}
  if ("Y" %in% Y$info[["chr"]]) {Y$info[Y$info[["chr"]] == "Y", "chr"] <- "24"}

  # Quarantee that there are no duplicated rows (probes) in the data
  # TODO: unless segmented data is used
  # (which should be explicitly indicated: add the option to function call)
  dupl <- duplicated(X$data)
  if (any(dupl)) {
    cat("Removing duplicate probe signals on X data..\n")
    X$data <- X$data[!dupl, ]
    X$info <- X$info[!dupl, ]
  }

  dupl <- duplicated(Y$data)
  if (any(dupl)) {
    cat("Removing duplicate probe signals on Y data..\n")
    Y$data <- Y$data[!dupl, ]
    Y$info <- Y$info[!dupl, ]
  }

  
  # First order chromosomes 1...22, X, Y, then chromosomes with other names
  if (is.null(chrs)) {
    chrs <- c(as.character(1:24), sort(setdiff(unique(X$info[["chr"]]), as.character(1:24))))
  } else {
    chrs <- as.character(chrs)
  }

  # If location information ('loc') for the probe is missing
  # but start and end positions are available, use them to calculate
  # probe middle bp location
  if (!"loc" %in% colnames(Y$info)) {
    # sometimes factors, sometimes numerics given; this should handle both cases correctly
    Y$info[["loc"]] <- (as.numeric(as.character(Y$info[, "start"])) + as.numeric(as.character(Y$info[, "end"])))/2
  }
  if (!"loc" %in% colnames(X$info)) {
    X$info[["loc"]] <- (as.numeric(as.character(X$info[, "start"])) + as.numeric(as.character(X$info[, "end"])))/2
  }
  
  message("Matching probes between the data sets..")
  xindices <- yindices <- vector()
  for (chr in chrs){
    # Note: chromosome arm information is used in the matching if it is available
    tmp <- get.neighboring.probes(X, Y, chr, max.dist)
    xindices <- c(xindices, tmp$xinds)
    yindices <- c(yindices, tmp$yinds)
  }

  xdat <- as.matrix(X$data[xindices,], length(xindices))
  ydat <- as.matrix(Y$data[yindices,], length(yindices))

  # Sometimes file reading (with csv at least) leads to situation where the last column is NA.
  # To avoid this and other cases, remove 'NA samples'.
  nainds <- (colMeans(is.na(xdat)) == 1 | colMeans(is.na(ydat)) == 1)
  if (sum(nainds) > 0) {
    xdat <- xdat[, !nainds]
    ydat <- ydat[, !nainds]
    warning(paste("Samples ", colnames(X$data)[nainds], " contained exclusively NA's; removed."))
  }
    
  newX <- list(data = xdat, info = X$info[xindices,])
  newY <- list(data = ydat, info = Y$info[yindices,])

  list(X = newX, Y = newY)

}

permute_dr_corr <- function(gepos, cnpos, method, ged, cnd){
  drcorrelate(matrix(ged[gepos,], nrow = 1),
              matrix(cnd[cnpos,], nrow = 1),
              method = method)
}


roc <- function (ordered.results, P) {

  # Calculate ROC curve
  
        #ordered results: best to worst
    #P: known positives
    #output: true positive rate and false positive rate
    
    #Check that all known positives are included in the original analysis i.e. ordered results list
    #if (!all(P %in% ordered.results)) {print("Warning: not all known positives are in the results list. Only included positives are used.")}
    positives<-P[P %in% ordered.results]    
    
    #Number of retrieved known cytobands
    N<-length(ordered.results) #total number of samples
    Np<-length(positives) #number of positives
    Nn<-N-Np #number of negatives

    TP<-cumsum(ordered.results %in% positives)
    FP<-cumsum(!(ordered.results %in% positives))
    tpr<-TP/Np #TP/(TP + FN) = TP.simCCA / P
    fpr<-FP/Nn #FP/(FP + TN) = FP.simCCA / N.simCCA

    list(tpr=tpr,fpr=fpr)
}


sym2gid <- function (gsym, xx, uniq = TRUE) {
  #symbol2geneid

  #symbol2geneid
  # Get mappings
  #if (is.null(xx)) {require("org.Hs.eg.db"); xx <- as.list(org.Hs.egALIAS2EG)} 
  # Remove NAs
  xx <- xx[!is.na(xx)]
  # pick just those in our data
  xx <- xx[names(xx) %in% gsym | names(xx) %in% tolower(gsym)]
	    	     
  gids <- list()
  for (sym in gsym) {
    #print(sym)
    if (sym %in% names(xx)) {
      gid <- unique(unlist(xx[sym]))
    } else if (tolower(sym) %in% names(xx)) {
      gid <- unique(unlist(xx[tolower(sym)]))
    } else {gid <- NULL}
    gids[[sym]] <- gid
  }
  # remove symbols that map to >1 geneids, 
  if (uniq) {
    gids <- gids[sapply(gids, length) == 1]
    return(gids)
  } else {return(gids)}
}

unitscale <- function(X, rm.na = TRUE, sd.value = NULL) {

  # scale each column to unit variance
  # remove col means from matrix X

  if (rm.na) {
    X2 <- matrix(NA,nrow=nrow(X),ncol=ncol(X))
    for (i in 1:ncol(X)) {
      x <- X[,i]
      nainds <- is.na(x)
      x <- x[!nainds]
      X2[!nainds,i] <- x/sd(x)
    }
  }
  if (!rm.na) {
    X2 <- apply(X,2,function(x){x/sd(x)})
  }
  if (length(sd.value)>0) {
    # Scale to predefined sd
    X2 <- apply(X2,2,function(x){x*sd.value})
  }

  dimnames(X2) <- dimnames(X)
  
  X2
}


intCNGEan.test.mod <- function (data.tuned, analysis.type, test.statistic, nperm = 10000, 
    eff.p.val.thres = 0.1) 
{
    uni.an.wcvm <- function(data.both, nosamp, a, nperm, low.ci.thres = 0.1) {
        wcvm.obs <- apply(data.both, 1, "wcvm.test.stats", nosamp, 
            a)
        data.perm <- data.both
        total.genes.on.chr <- dim(data.both)[1]
        remainders <- c(1:dim(data.perm)[1])
        steps <- sort(unique(c(0, 25, 50, 100, 150, 200, 250, 
            500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 
            2750, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 
            7000, 7500, 8000, 8500, 9000, 9500, seq(from = 10000, 
                to = 50000, by = 1000), nperm)))
        steps <- steps[steps <= nperm]
        wcvm.mat <- c()
        for (j in 1:(length(steps) - 1)) {
            for (i in 1:(steps[j + 1] - steps[j])) {
                if (((steps[j] + i)%%100) == 0) {
                  cat(paste(steps[j] + i, " of ", steps[length(steps)], 
                    " permutations done, and counting...", sep = ""), 
                    "\n")
                }
                x <- sample(1:nosamp, nosamp) + a * nosamp
                data.ran <- cbind(data.perm[, c(1:(a * nosamp))], 
                  data.perm[, x])
                wcvm.ran <- apply(data.ran, 1, "wcvm.test.stats", 
                  nosamp, a)
                wcvm.mat <- cbind(wcvm.mat, wcvm.ran)
            }
            perm.and.obs <- cbind(wcvm.mat, wcvm.obs)
            pvals <- apply(perm.and.obs, 1, "countth")/steps[j + 
                1]
            pbound <- sapply(pvals, "pvalbound", steps[j + 1])
            index <- cbind(1:length(pbound), pbound)
            pind <- index[pbound < low.ci.thres, 1]
            if (length(pind) > 2) {
                remainders <- remainders[pind]
            }
            else {
                pind <- c(1:length(remainders))
            }
            wcvm.mat <- wcvm.mat[pind, ]
            wcvm.obs <- wcvm.obs[pind]
            cat(paste(steps[j] + i, "of", steps[length(steps)], 
                " permutations done, and", length(remainders), 
                "of", total.genes.on.chr, "genes remaining...", 
                sep = " "), "\n")
            data.perm <- data.both[remainders, , drop = FALSE]
        }
        raw.pvals <- cbind(c(1:dim(data.both)[1]), rep(1, dim(data.both)[1]))
        raw.pvals[remainders, 2] <- pvals[pind]
        adjpvals <- cbind(raw.pvals, p.adjust(raw.pvals[, 2], 
            "BH"))
        adjpvals[, 2:3] <- round(adjpvals[, 2:3], digits = 4)
        colnames(adjpvals) <- c("clone.id", "raw.p", "adj.p")
        return(adjpvals)
    }
    uni.an.prob <- function(data.both, nosamp, a, nperm, low.ci.thres = 0.1) {
        prob.obs <- apply(data.both, 1, "prob.test.stats", nosamp, 
            a)
        data.perm <- data.both
        total.genes.on.chr <- dim(data.both)[1]
        remainders <- c(1:dim(data.perm)[1])
        steps <- sort(unique(c(0, 25, 50, 100, 150, 200, 250, 
            500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 
            2750, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 
            7000, 7500, 8000, 8500, 9000, 9500, seq(from = 10000, 
                to = 50000, by = 1000), nperm)))
        steps <- steps[steps <= nperm]
        prob.mat <- c()
        for (j in 1:(length(steps) - 1)) {
            for (i in 1:(steps[j + 1] - steps[j])) {
                if (((steps[j] + i)%%100) == 0) {
                  cat(paste(steps[j] + i, " of ", steps[length(steps)], 
                    " permutations done, and counting...", sep = ""), 
                    "\n")
                }
                x <- sample(1:nosamp, nosamp) + a * nosamp
                data.ran <- cbind(data.perm[, c(1:(a * nosamp))], 
                  data.perm[, x])
                prob.ran <- apply(data.ran, 1, "prob.test.stats", 
                  nosamp, a)
                prob.mat <- cbind(prob.mat, prob.ran)
            }
            perm.and.obs <- cbind(prob.mat, prob.obs)
            pvals <- apply(perm.and.obs, 1, "countth")/steps[j + 
                1]
            pbound <- sapply(pvals, "pvalbound", steps[j + 1])
            index <- cbind(1:length(pbound), pbound)
            pind <- index[pbound < low.ci.thres, 1]
            if (length(pind) > 2) {
                remainders <- remainders[pind]
            }
            else {
                pind <- c(1:length(remainders))
            }
            prob.mat <- prob.mat[pind, ]
            prob.obs <- prob.obs[pind]
            cat(paste(steps[j] + i, "of", steps[length(steps)], 
                " permutations done, and", length(remainders), 
                "of", total.genes.on.chr, "genes remaining...", 
                sep = " "), "\n")
            data.perm <- data.both[remainders, , drop = FALSE]
        }
        raw.pvals <- cbind(c(1:dim(data.both)[1]), rep(1, dim(data.both)[1]))
        raw.pvals[remainders, 2] <- pvals[pind]
        adjpvals <- cbind(raw.pvals, p.adjust(raw.pvals[, 2], 
            "BH"))
        adjpvals[, 2:3] <- round(adjpvals[, 2:3], digits = 4)
        colnames(adjpvals) <- c("clone.id", "raw.p", "adj.p")
        return(adjpvals)
    }
    R2.stat <- function(data.both, a, nosamp) {
        r2.numerator <- function(data.both, a, nosamp) {
            alpha.ind <- matrix(data.both[c(1:(a * nosamp))], 
                ncol = a, byrow = TRUE)
            alpha.1.mom <- apply(alpha.ind, 2, mean)
            alpha.2.mom <- t(alpha.ind) %*% alpha.ind/nosamp
            call.means <- NULL
            for (c in 1:a) {
                if (c == 1) {
                  c.help <- 2
                }
                else {
                  c.help <- 1
                }
                c.correction <- (alpha.2.mom[c, c]/alpha.2.mom[1, 
                  2] - alpha.1.mom[c]/alpha.1.mom[c.help])^(-1)
                slh <- as.numeric((data.both[a * (c(1:nosamp) - 
                  1) + c]/alpha.2.mom[1, 2] - 1/alpha.1.mom[c.help]) * 
                  data.both[a * nosamp + c(1:nosamp)])
                call.means <- cbind(call.means, c.correction * 
                  sum(slh)/nosamp)
            }
            return(sum((data.both[a * nosamp + c(1:nosamp)] - 
                alpha.ind %*% t(call.means))^2)/(nosamp - 1))
        }
        r2.denominator <- apply(data.both[, c((a * nosamp + 1):((a + 
            1) * nosamp))], 1, var)
        numerator <- apply(data.both, 1, "r2.numerator", a = a, 
            nosamp = nosamp)
        R2 <- 1 - numerator/r2.denominator
        for (i in 1:length(R2)) {
            R2[i] <- min(1, max(0, R2[i]))
        }
        return(R2)
    }
    shrin.an.prob <- function(data.both, nosamp, a, nperm, low.ci.thres = 0.1, 
        datacgh.org) {
        shrunken.prob.test.stats <- function(reg.bounds.lambda, 
            data.both, nosamp, a) {
            shrunken.test.stats.wrong.format <- apply(reg.bounds.lambda, 
                1, "shrunken.prob.test.stats.per.reg", cgh.em = data.both, 
                nosamp = nosamp, a = a)
            if (is.numeric(shrunken.test.stats.wrong.format)) {
                shrunken.test.stats <- shrunken.test.stats.wrong.format
            }
            if (is.matrix(shrunken.test.stats.wrong.format)) {
                shrunken.test.stats <- as.numeric(shrunken.test.stats.wrong.format)
            }
            if (is.list(shrunken.test.stats.wrong.format)) {
                shrunken.test.stats <- NULL
                for (i in 1:dim(reg.bounds.lambda)[1]) {
                  shrunken.test.stats <- c(shrunken.test.stats, 
                    shrunken.test.stats.wrong.format[[i]])
                }
            }
            return(shrunken.test.stats)
        }
        shrunken.prob.test.stats.per.reg <- function(bounds.lambda, 
            cgh.em, nosamp, a) {
            lambda <- bounds.lambda[3]
            bounds <- bounds.lambda[c(1:2)]
            if (bounds[1] != bounds[2]) {
                cgh.em <- cgh.em[c(bounds[1]:bounds[2]), ]
                marg.test.stats <- apply(cgh.em, 1, "prob.test.stats", 
                  nosamp, a)
                shrunken.test.stats <- lambda * marg.test.stats + 
                  (1 - lambda) * mean(marg.test.stats)
            }
            else {
                shrunken.test.stats <- wcvm.test.stats(cgh.em[bounds[1], 
                  ], nosamp, a)
            }
            return(shrunken.test.stats)
        }
        cat("construct regions...", "\n")
        reg.bounds <- find.cgh.reg.data(data.both, a, nosamp, 
            datacgh.org)
        cat("calculate shrinkage parameters...", "\n")
        lambda.of.reg <- apply(reg.bounds, 1, "lambda.per.reg", 
            data.both = data.both, a = a, nosamp = nosamp)
        reg.bounds.lambda <- cbind(reg.bounds, lambda.of.reg)
        colnames(reg.bounds.lambda) <- NULL
        cat("calculate observed test statistics...", "\n")
        shrunken.prob.obs <- shrunken.prob.test.stats(reg.bounds.lambda, 
            data.both, nosamp, a)
        cat("calculate null distribution...", "\n")
        prob.obs <- shrunken.prob.obs
        data.perm <- data.both
        total.genes.on.chr <- dim(data.both)[1]
        remainders <- c(1:dim(data.perm)[1])
        steps <- sort(unique(c(0, 25, 50, 100, 150, 200, 250, 
            500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 
            2750, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 
            7000, 7500, 8000, 8500, 9000, 9500, seq(from = 10000, 
                to = 50000, by = 1000), nperm)))
        steps <- steps[steps <= nperm]
        prob.mat <- c()
        for (j in 1:(length(steps) - 1)) {
            for (i in 1:(steps[j + 1] - steps[j])) {
                if (((steps[j] + i)%%100) == 0) {
                  cat(paste(steps[j] + i, " of ", steps[length(steps)], 
                    " permutations done, and counting...", sep = ""), 
                    "\n")
                }
                x <- sample(1:nosamp, nosamp) + a * nosamp
                data.ran <- cbind(data.perm[, c(1:(a * nosamp))], 
                  data.perm[, x])
                prob.ran <- shrunken.prob.test.stats(reg.bounds.lambda, 
                  data.ran, nosamp, a)
                prob.mat <- cbind(prob.mat, prob.ran[remainders])
            }
            perm.and.obs <- cbind(prob.mat, prob.obs)
            pvals <- apply(perm.and.obs, 1, "countth")/steps[j + 
                1]
            pbound <- sapply(pvals, "pvalbound", steps[j + 1])
            index <- cbind(1:length(pbound), pbound)
            pind <- index[pbound < low.ci.thres, 1]
            if (length(pind) > 2) {
                remainders <- remainders[pind]
            }
            else {
                pind <- c(1:length(remainders))
            }
            prob.mat <- prob.mat[pind, , drop = FALSE]
            prob.obs <- prob.obs[pind]
            cat(paste(steps[j] + i, "of", steps[length(steps)], 
                " permutations done, and", length(remainders), 
                "of", total.genes.on.chr, "genes remaining...", 
                sep = " "), "\n")
        }
        raw.pvals <- cbind(c(1:dim(data.both)[1]), rep(1, dim(data.both)[1]))
        raw.pvals[remainders, 2] <- pvals[pind]
        adjpvals <- cbind(raw.pvals, p.adjust(raw.pvals[, 2], 
            "BH"))
        adjpvals[, 2:3] <- round(raw.pvals[, 2:3], digits = 4)
        reg.details <- NULL
        for (i in 1:dim(reg.bounds.lambda)[1]) {
            reg.length <- (reg.bounds.lambda[i, 2] - reg.bounds.lambda[i, 
                1] + 1)
            reg.details <- rbind(reg.details, cbind(rep(i, reg.length), 
                matrix(rep(reg.bounds.lambda[i, ], reg.length), 
                  nrow = reg.length, byrow = TRUE)))
        }
        adjpvals <- cbind(adjpvals[, 1], reg.details, adjpvals[, 
            2:3])
        colnames(adjpvals) <- c("clone.id", "reg.id", "begin.reg", 
            "end.reg", "shrinkage", "raw.p", "adj.p")
        return(adjpvals)
    }
    shrin.an.wcvm <- function(data.both, nosamp, a, nperm, low.ci.thres = 0.1, 
        datacgh.org) {
        shrunken.wcvm.test.stats <- function(reg.bounds.lambda, 
            data.both, nosamp, a) {
            shrunken.test.stats.wrong.format <- apply(reg.bounds.lambda, 
                1, "shrunken.wcvm.test.stats.per.reg", cgh.em = data.both, 
                nosamp = nosamp, a = a)
            if (is.numeric(shrunken.test.stats.wrong.format)) {
                shrunken.test.stats <- shrunken.test.stats.wrong.format
            }
            if (is.matrix(shrunken.test.stats.wrong.format)) {
                shrunken.test.stats <- as.numeric(shrunken.test.stats.wrong.format)
            }
            if (is.list(shrunken.test.stats.wrong.format)) {
                shrunken.test.stats <- NULL
                for (i in 1:dim(reg.bounds.lambda)[1]) {
                  shrunken.test.stats <- c(shrunken.test.stats, 
                    shrunken.test.stats.wrong.format[[i]])
                }
            }
            return(shrunken.test.stats)
        }
        shrunken.wcvm.test.stats.per.reg <- function(bounds.lambda, 
            cgh.em, nosamp, a) {
            lambda <- bounds.lambda[3]
            bounds <- bounds.lambda[c(1:2)]
            if (bounds[1] != bounds[2]) {
                cgh.em <- cgh.em[c(bounds[1]:bounds[2]), ]
                marg.test.stats <- apply(cgh.em, 1, "wcvm.test.stats", 
                  nosamp, a)
                shrunken.test.stats <- lambda * marg.test.stats + 
                  (1 - lambda) * mean(marg.test.stats)
            }
            else {
                shrunken.test.stats <- wcvm.test.stats(cgh.em[bounds[1], 
                  ], nosamp, a)
            }
            return(shrunken.test.stats)
        }
        cat("construct regions...", "\n")
        reg.bounds <- find.cgh.reg.data(data.both, a, nosamp, 
            datacgh.org)
        cat("calculate shrinkage parameters...", "\n")
        lambda.of.reg <- apply(reg.bounds, 1, "lambda.per.reg", 
            data.both = data.both, a = a, nosamp = nosamp)
        reg.bounds.lambda <- cbind(reg.bounds, lambda.of.reg)
        colnames(reg.bounds.lambda) <- NULL
        cat("calculate observed test statistics...", "\n")
        shrunken.wcvm.obs <- shrunken.wcvm.test.stats(reg.bounds.lambda, 
            data.both, nosamp, a)
        cat("calculate null distribution...", "\n")
        wcvm.obs <- shrunken.wcvm.obs
        data.perm <- data.both
        total.genes.on.chr <- dim(data.both)[1]
        remainders <- c(1:dim(data.perm)[1])
        steps <- sort(unique(c(0, 25, 50, 100, 150, 200, 250, 
            500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 
            2750, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 
            7000, 7500, 8000, 8500, 9000, 9500, seq(from = 10000, 
                to = 50000, by = 1000), nperm)))
        steps <- steps[steps <= nperm]
        wcvm.mat <- c()
        for (j in 1:(length(steps) - 1)) {
            for (i in 1:(steps[j + 1] - steps[j])) {
                if (((steps[j] + i)%%100) == 0) {
                  cat(paste(steps[j] + i, " of ", steps[length(steps)], 
                    " permutations done, and counting...", sep = ""), 
                    "\n")
                }
                x <- sample(1:nosamp, nosamp) + a * nosamp
                data.ran <- cbind(data.perm[, c(1:(a * nosamp))], 
                  data.perm[, x])
                wcvm.ran <- shrunken.wcvm.test.stats(reg.bounds.lambda, 
                  data.ran, nosamp, a)
                wcvm.mat <- cbind(wcvm.mat, wcvm.ran[remainders])
            }
            perm.and.obs <- cbind(wcvm.mat, wcvm.obs)
            pvals <- apply(perm.and.obs, 1, "countth")/steps[j + 
                1]
            pbound <- sapply(pvals, "pvalbound", steps[j + 1])
            index <- cbind(1:length(pbound), pbound)
            pind <- index[pbound < low.ci.thres, 1]
            if (length(pind) > 2) {
                remainders <- remainders[pind]
            }
            else {
                pind <- c(1:length(remainders))
            }
            wcvm.mat <- wcvm.mat[pind, , drop = FALSE]
            wcvm.obs <- wcvm.obs[pind]
            cat(paste(steps[j] + i, "of", steps[length(steps)], 
                " permutations done, and", length(remainders), 
                "of", total.genes.on.chr, "genes remaining...", 
                sep = " "), "\n")
        }
        raw.pvals <- cbind(c(1:dim(data.both)[1]), rep(1, dim(data.both)[1]))
        raw.pvals[remainders, 2] <- pvals[pind]
        adjpvals <- cbind(raw.pvals, p.adjust(raw.pvals[, 2], 
            "BH"))
        adjpvals[, 2:3] <- round(raw.pvals[, 2:3], digits = 4)
        reg.details <- NULL
        for (i in 1:dim(reg.bounds.lambda)[1]) {
            reg.length <- (reg.bounds.lambda[i, 2] - reg.bounds.lambda[i, 
                1] + 1)
            reg.details <- rbind(reg.details, cbind(rep(i, reg.length), 
                matrix(rep(reg.bounds.lambda[i, ], reg.length), 
                  nrow = reg.length, byrow = TRUE)))
        }
        adjpvals <- cbind(adjpvals[, 1], reg.details, adjpvals[, 
            2:3])
        colnames(adjpvals) <- c("clone.id", "reg.id", "begin.reg", 
            "end.reg", "shrinkage", "raw.p", "adj.p")
        return(adjpvals)
    }
    lambda.per.reg <- function(bounds, data.both, a, nosamp) {
        if (bounds[1] != bounds[2]) {
            data.both.org <- matrix(data.both[c(bounds[1]:bounds[2]), 
                ], nrow = (bounds[2] - bounds[1] + 1))
            data.both.org <- data.both.org[, c((dim(data.both)[2] - 
                dim(data.both)[2]/3 + 1):dim(data.both)[2])]
            slh <- cor(t(data.both.org), method = "spearman")
            lambda.final <- 1 - max(0, mean(slh[upper.tri(slh)]))
        }
        if (bounds[1] == bounds[2]) {
            lambda.final <- 1
        }
        return(lambda.final)
    }
    find.cgh.reg.data <- function(data.both, a, nosamp, datacgh.org) {
        cgh.data <- data.both[, c(1:(a * nosamp))]
        splitter <- list()
        splitter[[1]] <- c(1)
        index.temp <- 1
        j <- 1
        for (i in 1:(dim(cgh.data)[1] - 1)) {
            if (all(cgh.data[i, ] == cgh.data[i + 1, ])) {
                index.temp <- c(index.temp, i + 1)
                splitter[[j]] <- index.temp
            }
            else {
                index.temp <- i + 1
                j <- j + 1
                splitter[[j]] <- index.temp
            }
        }
        region.details <- NULL
        for (i in 1:length(splitter)) {
            region.details <- rbind(region.details, c(min(splitter[[i]]), 
                max(splitter[[i]])))
        }
        return(region.details)
    }
    pvalbound <- function(pval, np) {
        return(pval - sqrt(pval * (1 - pval)/np) * 3.09)
    }
    countth.eff <- function(statlist, threshold) {
        return(length(statlist[statlist >= threshold]))
    }
    countth <- function(statlist) {
        threshold <- as.numeric(statlist[length(statlist)])
        statlist <- as.numeric(statlist[c(1:(length(statlist) - 
            1))])
        return(length(statlist[statlist >= threshold]))
    }
    pval.perm.marg <- function(observed, permuted, nperm) {
        perm.and.obs <- cbind(permuted, observed)
        return(apply(perm.and.obs, 1, "countth")/nperm)
    }
    rawps <- function(stats.obs, nulldists, nperm) {
        pval.ln <- pval.perm.marg(stats.obs, nulldists, nperm)
        return(pval.ln)
    }
    wcvm.test.stats <- function(cgh.em, nosamp, a) {
        cgh.2cat <- matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE)
        alphaas <- t(cgh.2cat) %*% cgh.2cat
        cs <- as.numeric(solve(alphaas) %*% matrix(c(-1, 1), 
            ncol = 1))
        cgh.em <- cbind(cgh.2cat, cgh.em[c((a * nosamp + 1):((a + 
            1) * nosamp))])
        cgh.em <- cbind(cgh.em[order(cgh.em[, 3]), ], rep(1/dim(cgh.em)[1], 
            dim(cgh.em)[1]))
        cgh.em <- cbind(cgh.em, cumsum(cgh.em[, 4]), cumsum(cgh.em[, 
            1] * cgh.em[, dim(cgh.em)[2]]) * cs[1], cumsum(cgh.em[, 
            2] * cgh.em[, dim(cgh.em)[2]]) * cs[2])
        test.stat <- -sum(cgh.em[, 7] + cgh.em[, 6])/nosamp
        return(test.stat)
    }
    wmw.test.stats <- function(cgh.em, nosamp, a) {
        cgh.em <- cbind(matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE), cgh.em[c((a * nosamp + 1):((a + 1) * 
            nosamp))])
        cgh.em <- cgh.em[order(cgh.em[, (a + 1)]), ]
        cgh.em <- cbind(cgh.em, rbind(rep(0, a), apply(cgh.em[, 
            1:a], 2, cumsum)[-nosamp, ]))
        test.stat <- sum(cgh.em[, a + 2] * cgh.em[, 2])
        return(test.stat)
    }
    prob.test.stats <- function(cgh.em, nosamp, a) {
        cgh.2cat <- matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE)
        alphaas <- t(cgh.2cat) %*% cgh.2cat
        cs <- c(det(alphaas), alphaas[1, 2] * sum(alphaas)/2)
        cgh.em <- cbind(matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE), cgh.em[c((a * nosamp + 1):((a + 1) * 
            nosamp))])
        cgh.em <- cgh.em[order(cgh.em[, (a + 1)]), ]
        cgh.em <- cbind(cgh.em, rbind(rep(0, a), apply(cgh.em[, 
            1:a], 2, cumsum)[-nosamp, ]))
        test.stat <- (sum(cgh.em[, a + 2] * cgh.em[, 2]) - cs[2])/cs[1]
        return(test.stat)
    }
    results <- NULL
    for (chr in 1:length(unique(data.tuned$ann[, 1]))) {
        chrs.present <- unique(data.tuned$ann[, 1])
        set.seed(7396)
        genes.per.time <- which(data.tuned$ann[, 1] == chrs.present[chr])
        data.both <- data.tuned$datafortest[genes.per.time, , 
            drop = FALSE]
        if (dim(data.both)[1] == 1) {
            no.genes.only.one <- TRUE
            data.both <- rbind(data.both, data.both)
        }
        else {
            no.genes.only.one <- FALSE
        }
        cat(paste("chromosome ", chrs.present[chr], " started...", 
            sep = ""), "\n")
        cat(paste(length(genes.per.time), " genes to be tested...", 
            sep = ""), "\n")
        if (analysis.type == "univariate") {
            if (test.statistic == "wcvm") {
                adjpvals <- uni.an.wcvm(data.both, data.tuned$nosamp, 
                  2, nperm, low.ci.thres = eff.p.val.thres)
            }
            if (test.statistic == "wmw") {
                adjpvals <- uni.an.prob(data.both, data.tuned$nosamp, 
                  2, nperm, low.ci.thres = eff.p.val.thres)
            }
            results.temp <- data.frame()
            R2 <- R2.stat(data.both, 2, data.tuned$nosamp)
            R2[is.na(data.tuned$alleffects[genes.per.time])] <- NA
            if (no.genes.only.one) {
                results.temp <- c(data.tuned$genestotest[genes.per.time], 
                  data.tuned$lossorgain[genes.per.time], round(data.tuned$callprobs[genes.per.time, 
                    ], digits = 4), round(data.tuned$alleffects[genes.per.time], 
                    digits = 4), round(R2[1], digits = 4), round(adjpvals[1, 
                    2:3], digits = 4))
                names(results.temp)[1:6] <- c("gene.id", "comparison", 
                  "av.probs.1", "av.probs.2", "effect.size", 
                  "R2")
                results <- rbind(results, c(data.tuned$ann[genes.per.time, 
                  ], results.temp))
            }
            else {
                results.temp <- cbind(data.tuned$genestotest[genes.per.time], 
                  data.tuned$lossorgain[genes.per.time], round(data.tuned$callprobs[genes.per.time, 
                    ], digits = 4), round(data.tuned$alleffects[genes.per.time], 
                    digits = 4), round(R2, digits = 4), round(adjpvals[, 
                    2:3], digits = 4))
                colnames(results.temp)[1:6] <- c("gene.id", "comparison", 
                  "av.probs.1", "av.probs.2", "effect.size", 
                  "R2")
                results <- rbind(results, cbind(data.tuned$ann[genes.per.time, 
                  ], results.temp))
            }
        }
        if (analysis.type == "regional") {
            if (test.statistic == "wcvm") {
                adjpvals <- shrin.an.wcvm(data.both, data.tuned$nosamp, 
                  2, nperm, low.ci.thres = eff.p.val.thres, datacgh.org = data.both[, 
                    c(1:(2 * data.tuned$nosamp))])
            }
            if (test.statistic == "wmw") {
                adjpvals <- shrin.an.prob(data.both, data.tuned$nosamp, 
                  2, nperm, low.ci.thres = eff.p.val.thres, datacgh.org = data.both[, 
                    c(1:(2 * data.tuned$nosamp))])
            }
            results.temp <- data.frame()
            R2 <- R2.stat(data.both, 2, data.tuned$nosamp)
            R2[is.na(data.tuned$alleffects[genes.per.time])] <- NA
            if (no.genes.only.one) {
                results.temp <- c(data.tuned$genestotest[genes.per.time], 
                  data.tuned$lossorgain[genes.per.time], round(data.tuned$callprobs[genes.per.time, 
                    ], digits = 4), round(data.tuned$alleffects[genes.per.time], 
                    digits = 4), round(R2[1], digits = 4), adjpvals[1, 
                    -1])
                names(results.temp)[1:6] <- c("gene.id", "comparison", 
                  "av.probs.1", "av.probs.2", "effect.size", 
                  "R2")
                results.temp[12] <- round(results.temp[12], digits = 4)
                results.temp[11] <- round(results.temp[11], digits = 4)
                results.temp[9] <- results.temp[8]
                results <- rbind(results, c(data.tuned$ann[genes.per.time, 
                  ], results.temp))
            }
            else {
                results.temp <- cbind(data.tuned$genestotest[genes.per.time], 
                  data.tuned$lossorgain[genes.per.time], round(data.tuned$callprobs[genes.per.time, 
                    ], digits = 4), round(data.tuned$alleffects[genes.per.time], 
                    digits = 4), round(R2, digits = 4), adjpvals[, 
                    -1])
                colnames(results.temp)[1:6] <- c("gene.id", "comparison", 
                  "av.probs.1", "av.probs.2", "effect.size", 
                  "R2")
                results.temp[, 12] <- round(results.temp[, 12], 
                  digits = 4)
                results.temp[, 11] <- round(results.temp[, 11], 
                  digits = 4)
                results <- rbind(results, cbind(data.tuned$ann[genes.per.time, 
                  ], results.temp))
            }
        }
        cat(paste("chromosome ", chrs.present[chr], " done...", 
            sep = ""), "\n")
    }
    results[, dim(results)[2]] <- p.adjust(results[, dim(results)[2] - 
        1], "none")
    cat("ready: testing done", "\n")
    return(results)
}
