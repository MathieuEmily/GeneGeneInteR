GBIGM.test <- function(Y, G1, G2, n.perm = 1000){

	Y.arg <- deparse(substitute(Y))
	G1.arg <- deparse(substitute(G1))
	G2.arg <- deparse(substitute(G2))
  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if(class(n.perm)!="numeric"){
    stop("n.boot must be numeric.")
  } else if(n.perm<1){
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

  n.ind<-nrow(X1)

  # R a une precision jusqu'a 2*3^33+1, donc au-dessus de 33 SNPs totaux sur les 2 genes,
  # il faut utiliser la fonction paste().

  if(ncol(X1)+ncol(X2)<34){
    ind34<-1
    G1<-base3to10(X1)
    G2<-base3to10(X2)
    G1.2<-apply(cbind(G1,3^ncol(X1)*G2),1,sum) # g?notype de l'individu X1-X2 transcrit en base10
  }else{
    ind34<-0
    G1<-apply(X1,1,paste,collapse=" ")
    G2<-apply(X2,1,paste,collapse=" ")
    G1.2<-apply(cbind(G1,G2),1,paste,collapse=" ")
  }

  HG1<-entropy.vec(G1)
  HG2<-entropy.vec(G2)
  HG1.2<-entropy.vec(G1.2)

  Delta1.2init<-InfoGainRat(Y,G1,G2,G1.2,HG1,HG2,HG1.2,ind34)
  D<-list()
  for(i in 1:n.perm){
    D[[i]]<-sample(Y,n.ind)
  }
  Delta1.2<-rep(NA,n.perm)
  Delta1.2<-lapply(1:n.perm,function(x){InfoGainRat(D[[x]],G1,G2,G1.2,HG1,HG2,HG1.2,ind34)})
  pval<-sum(Delta1.2>=Delta1.2init)/n.perm
	 stat <- Delta1.2init
	names(stat)="DeltaR1,2"
	# list.param <- list(n.perm=n.perm)
	# res <- list(statistic=stat,p.value=pval,method="Gene-based Information Gain Method",parameter=list.param)
	# class(res) <- "GGItest"
  # return(res)
#  return(pval)

    null.value <- NULL
#names(null.value) <- "DeltaR1,2"
estimate <- c(stat)
names(estimate) <- c("DeltaR1,2")
parameters <- n.perm
names(parameters) <- "n.perm"
	res <- list(
		null.value=null.value,
		alternative="two.sided",
		method="Gene-based interaction based on Gene-based Information Gain Method",
		estimate= estimate,
		data.name=paste(Y.arg," and  (",G1.arg," , ",G2.arg,")",sep=""),
		statistic=stat,
		p.value=pval,
		parameters=parameters)
	class(res) <- "htest"
  return(res)

}

base3to10 <- function(gene){
  n.SNP<-ncol(gene)
  long<-nrow(gene)
  tab<-matrix(NA,ncol=n.SNP,nrow=long)
  for (i in 1:n.SNP){
    tab[,i]<-3^(i-1)*gene[,i]
  }
  b10<-apply(tab,1,sum)
  return(b10)
}

entropy.vec <- function(X){
  pourc<-percentofX(X)
  H<-(-sum(pourc*log(pourc)))
  return(H)
}

percentofX <- function(X){
  n.ind<-length(X)
  b10<-as.factor(X)
  vect<-tabulate(b10)
  pourc<-vect/n.ind
  return(pourc)
}

InfoGainRat <- function(D,tb1,tb2,tb1.2,Hc1,Hc2,Hc1c2,ind34){
  if(!ind34){
    tabD1<-apply(cbind(D,tb1),1,paste,collapse=" ")
    tabD2<-apply(cbind(D,tb2),1,paste,collapse=" ")
    tabD1.2<-apply(cbind(D,tb1.2),1,paste,collapse=" ")
  }else{
    tabD1<-apply(cbind(D,2*tb1),1,sum)
    tabD2<-apply(cbind(D,2*tb2),1,sum)
    tabD1.2<-apply(cbind(D,2*tb1.2),1,sum)
  }
  H1<-entropy.vec(tabD1)-Hc1
  H2<-entropy.vec(tabD2)-Hc2
  H1.2<-entropy.vec(tabD1.2)-Hc1c2
  if(min(H1,H2)==0){
    return(min(H1,H2)-H1.2)
  }else{
    return((min(H1,H2)-H1.2)/min(H1,H2))
  }
}
