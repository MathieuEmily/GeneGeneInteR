CCA.test <- function(Y, G1, G2, n.boot = 500){

	Y.arg <- deparse(substitute(Y))
	G1.arg <- deparse(substitute(G1))
	G2.arg <- deparse(substitute(G2))
  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if(class(n.boot)!="numeric"){
    stop("n.boot must be numeric.")
  } else if(n.boot<1){
    stop("n.boot must be strictly superior to 0.")
  } else if(nlevels(as.factor(Y))!=2){
    stop("Y must be a factor with 2 levels, most likely 0 and 1.")
  } else if(class(G1)!="SnpMatrix"){
    stop("G1 must be a SnpMatrix object.")
  } else if(class(G2)!="SnpMatrix"){
    stop("G2 must be a SnpMatrix object")
  } else if(nrow(G1)!=nrow(G2)){
    stop("Both G1 and G2 must contain the same number of individuals.")
  } else if(length(Y)!=nrow(G1)){
    stop("Y and both SnpMatrix objects must contain the same number of individuals.")
  } else if (sum(is.na(G1))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(G2))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(Y))!=0) {
    stop("The response variable must be complete. No NAs are allowed.")
  }

  X1 <- as(G1, "numeric")
  X2 <- as(G2, "numeric")

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  result<-"empty"
  try(result <- get.U(Y=Y,X1=X1,X2=X2,n.boot=n.boot))

  if(is.na(result[1])){
  	stat <- NA
    pval<-1

  } else {
    if(result[1]=="empty"){stop("P-value could not be computed, test statistic is missing")}
    stat <- result$stat
    pval <- 2*(1-pnorm(abs(stat)))
	z0 <- result$z0
	z1 <- result$z1
  }
	names(stat)="CCU"

	# list.param<-list(n.boot = n.boot)
	# res <- list(statistic=stat,p.value=pval,method="Canonical Correlation Analysis",parameter=list.param)
	# class(res) <- "GGItest"
  # return(res)
  
null.value <- 0
names(null.value) <- "CCU"
estimate <- c(z0,z1)
names(estimate) <- c("z0","z1")
parameters <- n.boot
names(parameters) <- "n.boot"
	res <- list(
		null.value=null.value,
		alternative="two.sided",
		method="Gene-based interaction based on Canonical Correspondance Analysis",
		estimate= estimate,
		data.name=paste(Y.arg," and  (",G1.arg," , ",G2.arg,")",sep=""),
		statistic=stat,
		p.value=pval,
		parameters=parameters)
	class(res) <- "htest"
  return(res)
#  return(comp.res$'Pr(>Chi)'[2])

  
}


get.U <- function(Y,X1,X2,n.boot=500){
  w0 <- which(Y==0)
  w1 <- which(Y==1)
  X1.0 <- as.matrix(X1[w0,])
  X2.0 <- as.matrix(X2[w0,])
  X1.1 <- as.matrix(X1[w1,])
  X2.1 <- as.matrix(X2[w1,])
  z0 <- get.z(X1.0,X2.0)
  z1 <- get.z(X1.1,X2.1)
  if(is.na(z0)&&is.na(z1)){
    warning("Canonical correlations between the first gene cases and the second gene cases equals to 1")
    warning("Canonical correlations between the first gene controls and the second gene controls equals to 1")
    return(0)
  }else{
    if(is.na(z0)){
      warning("Canonical correlations between the first gene controls and the second gene controls equals to 1")
      return(NA)
    }else{
      if(is.na(z1)){
        warning("Canonical correlations between the first gene cases and the second gene cases equals to 1")
        return(NA)
      }else{
        vz0<-NA;vz1<-NA;
        try(vz0 <- estim.var.z(X1.0,X2.0,n.boot=n.boot))
        try(vz1 <- estim.var.z(X1.1,X2.1,n.boot=n.boot))
        if(is.na(vz0)||is.na(vz1)){stop("The test statistic could not be computed, the variance estimator is missing")}
        else{return(list(stat=(z0-z1)/sqrt(vz0+vz1),z0=z0,z1=z1))}
      }
    }
  }
}


get.z <- function(X1,X2){
  tmp <- cancor(X1,X2)
  r <- tmp$cor[1]
  if(r>1-10^-16){
    z=NA;
  }else{
    z <- (1/2)*(log(1+r)-log(1-r))
  }
  return(z)
}

estim.var.z <- function(X1,X2,n.boot=500){
  z.vec <- rep(NA,times=n.boot)
  for (i in 1:n.boot){
    restart<-TRUE
    while(restart){
      ind.boot <- sample(1:nrow(X1),nrow(X1),replace=TRUE)
      X1.boot <- as.matrix(X1[ind.boot,])
      X2.boot <- as.matrix(X2[ind.boot,])
      z.vec[i] <- get.z(X1.boot,X2.boot)
      if(!is.na(z.vec[i])){restart<-FALSE}
    }
  }
  return(var(z.vec))
}
