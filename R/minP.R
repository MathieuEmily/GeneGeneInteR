minP.test <- function(Y, G1, G2){
	Y.arg <- deparse(substitute(Y))
	G1.arg <- deparse(substitute(G1))
	G2.arg <- deparse(substitute(G2))

  # Checking if gene splitting is needed
  if (ncol(G1) * ncol(G2) > 1000){
  	if (ncol(G1) > 30){
	  	## Searching for clusters in G1
  		distance.G1 <- snpStats::ld(G1, G1, stats='R.squared')
	  	distance.G1 <- as.dist(1 - distance.G1)
	  	clust.tree.G1 <- rioja::chclust(distance.G1)
	  	k1 <- cutree(clust.tree.G1, k=1:(ncol(G1)-30))
	  	max.G1 <- sapply(1:(ncol(G1)-30),FUN=function(i){return(max(table(as.factor(k1[,i]))))})
#	  	max.G1 <- max.G1[c(which(max.G1 > 20),which(max.G1 < 30)[1])]
	  	id.max.G1 <- which(max.G1 <= 30)[1]
 
  	} else{
  		k1  <- matrix(rep(1,ncol(G1)),ncol=1)
  		row.names(k1) <- colnames(G1)
  		max.G1 <- ncol(G1)
  		id.max.G1 <- 1
  		}
  	
  	
  	## Searching for clusters in G1
  	if (ncol(G2) >30){
  		distance.G2 <- snpStats::ld(G2, G2, stats='R.squared')
  		distance.G2 <- as.dist(1 - distance.G2)
  		clust.tree.G2 <- rioja::chclust(distance.G2)
  		k2 <- cutree(clust.tree.G2, k=1:(ncol(G2)-30))
  		max.G2 <- sapply(1:(ncol(G2)-30),FUN=function(i){return(max(table(as.factor(k2[,i]))))})
#  		max.G2 <- max.G2[c(which(max.G2 > 20),which(max.G2 < 30)[1])]
  		id.max.G2 <- which(max.G2 <= 30)[1]
	} else {
		k2  <- matrix(rep(1,ncol(G2)),ncol=1)
  		row.names(k2) <- colnames(G2)
  		max.G2 <- ncol(G2)
  		id.max.G2 <- 1
	}
	
	## Detection of the subset of SNPs that cut each gene in sufficiently small pieces 
	#mult <- outer(max.G1[ind.max.G1],max.G2[ind.max.G2],"*")
	#tmp.mult <- (mult < 900)*mult
	division <- c(id.max.G1,id.max.G2)
#	division <- which(tmp.mult==max(tmp.mult),arr.ind=TRUE)[1,]
#	division <- which(mult < 900,arr.ind=TRUE)[1,]



	boundaries.start.G1 <- c(1,1+as.numeric(which(sapply(1:(length(k1[,division[1]])-1),FUN=function(i){return(k1[,division[1]][i]!=k1[,division[1]][i+1])}))))
	boundaries.end.G1 <- c(as.numeric(which(sapply(1:(length(k1[,division[1]])-1),FUN=function(i){return(k1[,division[1]][i]!=k1[,division[1]][i+1])}))),ncol(G1))

	boundaries.start.G2 <- c(1,1+as.numeric(which(sapply(1:(length(k2[,division[2]])-1),FUN=function(i){return(k2[,division[2]][i]!=k2[,division[2]][i+1])}))))
	boundaries.end.G2 <- c(as.numeric(which(sapply(1:(length(k2[,division[2]])-1),FUN=function(i){return(k2[,division[2]][i]!=k2[,division[2]][i+1])}))),ncol(G2))
	} else{
		boundaries.start.G1 <- c(1)
		boundaries.start.G2 <- c(1)
		boundaries.end.G1 <- c(ncol(G1))
		boundaries.end.G2 <- c(ncol(G2))
	}
	
	sub.pairs.start <- expand.grid(boundaries.start.G1,boundaries.start.G2)
	sub.pairs.end <- expand.grid(boundaries.end.G1,boundaries.end.G2)
	sub.pairs <- data.frame(
		start.G1=sub.pairs.start[,1],
		end.G1=sub.pairs.end[,1],
		start.G2=sub.pairs.start[,2],
		end.G2=sub.pairs.end[,2]
	)
	
  # All pairs are iterated over
  pairs.p.val <- rep(NA,times=nrow(sub.pairs))
  pairs.stat <- rep(NA,times=nrow(sub.pairs))
  for (i in 1:(nrow(sub.pairs))){
    # Cluster boundaries for both genes
	G1.boundary <- (sub.pairs$start.G1[i]):(sub.pairs$end.G1[i])
	G2.boundary <- (sub.pairs$start.G2[i]):(sub.pairs$end.G2[i])
#    c.test <- list(p.value=runif(1),statistic=runif(1))
	c.test <-minP.test.2pairs(Y, G1[, G1.boundary], G2[, G2.boundary])
    pairs.p.val[i] <- c.test$p.value
    pairs.stat[i] <- c.test$statistic
    
  }


	tmp <- p.adjust(pairs.p.val, "BH")
	pval <- min(tmp)[1]
	stat <- pairs.stat[which.min(tmp)]
	names(stat)="Wmax"
#	res <- list(statistic=stat,p.value=pval,method="minP")
#	class(res) <- "GGItest"
 # return(res)
  
  	null.value <- NULL
#	names(null.value) <- "Wmax"
	estimate <- c(stat)
	names(estimate) <- c("Wmax")
	parameters <- NULL
#	names(parameters) <- ""
	res <- list(
		null.value=null.value,
		alternative="greater",
		method="Gene-based interaction based on minP method",
		estimate= estimate,
		data.name=paste(Y.arg," and  (",G1.arg," , ",G2.arg,")",sep=""),
		statistic=stat,
		p.value=pval,
		parameters=parameters)
	class(res) <- "htest"
  return(res)


}

minP.test.2pairs <- function(Y, G1, G2){
	Y.arg <- deparse(substitute(Y))
	G1.arg <- deparse(substitute(G1))
	G2.arg <- deparse(substitute(G2))

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if (nlevels(as.factor(Y)) != 2) {
    stop("response variable should be binary. (2 modes).")
  } else if (class(G1) != "SnpMatrix" | class(G2) != "SnpMatrix") {
    stop("G1 and G2 arguments should be SnpMatrix objects.")
  } else if (nrow(G1) != nrow(G2)) {
    stop("G1 and G2 should have same rows count.")
  } else if (length(Y) != nrow(G1) | length(Y) != nrow(G2)) {
    stop("Response variable should be conformant with genes matrices rows number.")
  } else if (sum(is.na(G1))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(G2))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(Y)) != 0) {
    stop("The response variable vector must be complete. No NAs are allowed.")
  }

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  # minP test is set up
  # SNP interactions are computed
  SSI.res <- SSI.test(Y, G1, G2)

  ## Computation of the correlation matrix
	MatCor1 <- 	snpStats::ld(G1, G1, stats="R")
	MatCor2 <- 	snpStats::ld(G2, G2, stats="R")
	n1 <- ncol(MatCor1)
	n2 <- ncol(MatCor2)
	n.pairs <- n1*n2
	sigma.matrix <- matrix(NA,ncol=n.pairs,nrow=n.pairs)
	for (i in 1:(n.pairs-1)){
		i1 <- floor((i-1)/n2)+1
		j1 <- i-(i1-1)*n2
		for (j in (i+1):n.pairs){
			i2 <- floor((j-1)/n2)+1
			j2 <- j-(i2-1)*n2
			sigma.matrix[i,j] <- sigma.matrix[j,i] <- MatCor1[i1,i2]*MatCor2[j1,j2]
		}
	}
	diag(sigma.matrix) <- 1

  # minP test is computed past this point
  SSI.min <- min(SSI.res, na.rm=TRUE)

  # Case in which 1 - SSI.min/2 == 1 (qnorm(1) = -Inf)
  if (SSI.min == 0) {
    GG.Pmin <- 0
    # Case in which 1 - SSI.min/2 == 1 (qnorm(1) = Inf)
  } else {
    lower <- rep(qnorm(SSI.min / 2), ncol(sigma.matrix))
    upper <- rep(qnorm(1 - SSI.min / 2), ncol(sigma.matrix))

    GG.Pmin <- 1 - mvtnorm::pmvnorm(lower=lower,upper=upper,sigma=sigma.matrix,maxpts=2.5E6,abseps = 0.01
    #abseps=1E-13
    )
  }
  if (attributes(GG.Pmin)$msg=="Completion with error > abseps"){
		warning(paste("p-values are approximated with error=",attributes(GG.Pmin)$error))
	}
	pval <- as.numeric(GG.Pmin)
	stat <- SSI.min
	names(stat)="Wmax"
	#res <- list(statistic=stat,p.value=pval,method="minP")
	#class(res) <- "GGItest"
#  return(res)

	null.value <- 0
	names(null.value) <- "Wmax"
	estimate <- c(stat)
	names(estimate) <- c("Wmax")
	parameters <- NULL
#	names(parameters) <- ""
	res <- list(
		null.value=null.value,
		alternative="greater",
		method="Gene-based interaction based on minP method",
		estimate= estimate,
		data.name=paste("Interaction between",G1.arg,"and",G2.arg,"in association with",Y.arg),
		statistic=stat,
		p.value=pval,
		parameters=parameters)
	class(res) <- "htest"
  return(res)


}
