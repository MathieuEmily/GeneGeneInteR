PLSPM.test <- function(Y, G1, G2, n.perm=500){

	Y.arg <- deparse(substitute(Y))
	G1.arg <- deparse(substitute(G1))
	G2.arg <- deparse(substitute(G2))
	

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if(nlevels(as.factor(Y))!=2){
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

  X <- cbind(X1,X2)

  Gene1 <- c(0,0)
  Gene2 <- c(1,0)

  w1 <- which(Y==1)
  w0 <- which(Y==0)

  XCases <- X[w1,]
  XControls <- X[w0,]

  my.path <- rbind(Gene1,Gene2)
  my.blocks <- list(1:ncol(X1),(ncol(X1)+1):(ncol(X1)+ncol(X2)))
  my.modes = c("A", "A")

  mod1<-NULL;
  mod0<-NULL;

  try(mod1 <- plspm::plspm(XCases,my.path,my.blocks, modes = my.modes), silent=TRUE)
  if(is.null(mod1)){
  	warning("P-value could not be computed. NA returned")
	list.param <- list(n.perm=n.perm)
  	res <- list(statistic=NA,p.value=NA,method="Partial Least Squares Path Modeling",parameter=list.param)
	class(res) <- "GGItest"
	return(res)
	}

  try(mod0 <- plspm::plspm(XControls,my.path,my.blocks, modes = my.modes),silent=TRUE)
  if(is.null(mod0)){
  	warning("P-value could not be computed. NA returned")
	list.param <- list(n.perm=n.perm)
  	res <- list(statistic=NA,p.value=NA,method="Partial Least Squares Path Modeling",parameter=list.param)
	class(res) <- "GGItest"
	return(res)
	}

  beta1 <- mod1$inner_model[[1]][2,1]
  vbeta1 <- mod1$inner_model[[1]][2,2]^2
  beta0 <- mod0$inner_model[[1]][2,1]
  vbeta0 <- mod0$inner_model[[1]][2,2]^2

  U <- (beta0-beta1)/sqrt(vbeta0+vbeta1)
  
  U.perm <- rep(NA,times=n.perm)
  for (i in 1:n.perm){
    restart<-TRUE
    while(restart){
	  	Y.perm <- sample(Y)
  		w1 <- which(Y.perm==1)
	  	w0 <- which(Y.perm==0)
  		XCases <- X[w1,]
	  	XControls <- X[w0,]
  		mod1<-NULL;
		mod0<-NULL;
		try(mod1 <- plspm::plspm(XCases,my.path,my.blocks, modes = my.modes), silent=TRUE)
#		if(is.null(mod1)){warning("P-value could not be computed. NA returned");return(NA)}
		try(mod0 <- plspm::plspm(XControls,my.path,my.blocks, modes = my.modes),silent=TRUE)
#		if(is.null(mod0)){warning("P-value could not be computed. NA returned");return(NA)}
		if (!is.null(mod1) & !is.null(mod0)){restart <- FALSE}
	}
	beta1 <- mod1$inner_model[[1]][2,1]
	vbeta1 <- mod1$inner_model[[1]][2,2]^2
	beta0 <- mod0$inner_model[[1]][2,1]
	vbeta0 <- mod0$inner_model[[1]][2,2]^2
	
	U.perm[i] <- (beta0-beta1)/sqrt(vbeta0+vbeta1)
  }
  
  #pval <- 2*(1-pnorm(abs(U)))
	  pval <- mean(abs(U.perm) > abs(U))
	  stat <- U
	names(stat)="U"
#	list.param <- list(n.perm=n.perm)
#	res <- list(statistic=stat,p.value=pval,method="Partial Least Squares Path Modeling",parameter=list.param)
#	class(res) <- "GGItest"
 # return(res)
  
    null.value <- 0
names(null.value) <- "U"
estimate <- c(beta0, beta1)
names(estimate) <- c("beta0","beta1")
parameters <- n.perm
names(parameters) <- "n.perm"
	res <- list(
		null.value=null.value,
		alternative="two.sided",
		method="Gene-based interaction based on Partial Least Squares Path Modeling",
		estimate= estimate,
		data.name=paste(Y.arg," and  (",G1.arg," , ",G2.arg,")",sep=""),
		statistic=stat,
		p.value=pval,
		parameters=parameters)
	class(res) <- "htest"
  return(res)


#  return(pval)
}
