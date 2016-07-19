tProd.test <- function(Y, G1, G2, tau = 0.05, n.sim = 1000){

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
  } else if(!is.numeric(tau) || tau < 0 || tau > 1){
    stop("tau must be a numeric in [0, 1].")
  } else if(!is.numeric(n.sim) || n.sim < 0){
    stop("n.sim must be a numeric strictly greater than 0.")
  }

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  # tTS test is set up
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

 ## Correlation is sorted according to the SSI tests.
	sigma.matrix[order(as.numeric(t(SSI.res))),order(as.numeric(t(SSI.res)))]


  # tTS test is computed past this point
  T0 <- tProd.stat(SSI.res,tau)

  sim <- mvtnorm::rmvnorm(n=n.sim, sigma = sigma.matrix, method="svd")
  Pvalsim <- 2*(1-pnorm(abs(sim)))
  Ti <- apply(Pvalsim,1,tProd.stat,tau)

	pval <- mean(Ti < T0,na.rm=TRUE)
	
	stat <- T0
	names(stat)="tProd"
	res <- list(statistic=stat,p.value=pval,method="Truncated Product")
	class(res) <- "GGItest"
  return(res)

#  return(sum(Ti[Ti<T0],na.rm=T)/n.sim)
}

# Function that computes test statistic for
# tProd method.
tProd.stat <- function(Pval, tau){
  Pval <- sort(as.vector(Pval))
  return(prod(Pval[Pval<tau]))
#  if(any(Pval<tau)){return(prod(Pval[Pval<tau]))}
#  else(NA)
}
