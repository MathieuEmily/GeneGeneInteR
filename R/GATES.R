gates.test <- function(Y, G1, G2, alpha = 0.05, me.est = c("ChevNy", "Keff", "LiJi","Galwey")){
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
  } else if (!is.numeric(alpha) | !(alpha >= 0 & alpha <= 1)) {
    stop("alpha argument should a numeric in [0, 1].")
  } else if (!is.character(me.est)) {
    stop("me.est argument should be a character string either ChevNy, Keff, LiJi or Galwey")
  }

  me.est <- try(match.arg(me.est))

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  if (class(me.est) == "try-error") {
    stop("me.est argument has an incorrect value. Select one: ChevNy, Keff, LiJi, Galwey.")
  }

  # GATES test is set up
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


  # GATES test is computed past this point
  sorted.SSI <- sort(SSI.res)

  me <- switch(me.est,
               ChevNy = ChevNy.me(sigma.matrix),
               Keff   = Keff.me(sigma.matrix, alpha),
               LiJi   = LiJi.me(sigma.matrix),
               Galwey   = Galwey.me(sigma.matrix)
  )

  GG.PGates <- min(me$Meff * sorted.SSI/me$mej)


	pval <- as.numeric(GG.PGates)
	stat <- sorted.SSI[which.min(me$Meff * sorted.SSI/me$mej)]
	names(stat)="GATES"
	res <- list(statistic=stat,p.value=pval,method="GATES")
	class(res) <- "GGItest"
  return(res)


#  return(GG.PGates)
}

# Cheverud-Nyholt Me estimation method
ChevNy.me <- function(sigma.matrix) {
  if (class(sigma.matrix) != "matrix") {
    stop("sigma.matrix argument should be a numeric matrix.")
  }

  N <- ncol(sigma.matrix)
  Meff <- 1 +(1/N)*sum(1-sigma.matrix^2, na.rm=TRUE)

  mej <- numeric(ncol(sigma.matrix))
  for (i in 1:ncol(sigma.matrix)){
    mej[i] <- 1 + (1/i)*sum(1-sigma.matrix[1:i, 1:i]^2, na.rm=TRUE)
  }

  return(list(Meff=Meff, mej=mej))
}

# Meff and Mej are estimated through the computation of
# Keff and kj
Keff.me <- function(sigma.matrix, alpha=0.05) {
  if (class(sigma.matrix) != "matrix") {
    stop("sigma.matrix argument should be a numeric matrix.")
  } else if (alpha < 0 | alpha > 1) {
    stop("alpha argument should be comprised between 0 and 1.")
  }

  N  <- ncol(sigma.matrix)
  kj <- rep(NA,N)
  kj[1] <- 0
  for (i in 2:length(kj)){
    if (alpha > 0.01) {
      rj    <- max(abs(sigma.matrix[(1:(i-1)), i]))
      sig   <- qnorm(1-(alpha/2))

      f     <- function(x, alpha, rj, sd) { exp(-0.5*x^2) * pnorm((rj*x - sd)/sqrt(1 - rj^2)) }
      int.f <- integrate(f, lower=-sig, upper=sig, alpha=alpha, rj=rj, sd=sig)$value

      kj[i] <- ( (1/log(1-alpha)) * log(1 - (1/(1 - alpha)) * sqrt(2/pi) * int.f) )
    } else {
      rj <- max(abs(sigma.matrix[(1:(i-1)), i]))
      kj[i] <- sqrt(1 - rj^(-1.31*log10(alpha)))
    }
  }

  Keff <- 1 + sum(kj)
  mej <- 1 + cumsum(kj)

  return(list(Meff=Keff, mej=mej))
}

# Li & Ji method to estimate Meff & mej
LiJi.me <- function(sigma.matrix) {
  if (class(sigma.matrix) != "matrix") {
    stop("sigma.matrix argument should be a numeric matrix.")
  }

  N <- ncol(sigma.matrix)

  f <- function(x) { ifelse(x >= 1, 1, 0) + (x - floor(x)) }

  res.vec=rep(NA,times=N)
  res.vec[1]=1
  for (i in 2:N){
  	eig.val <- eigen(sigma.matrix[1:i,1:i])$values
  	eig.val <- f(abs(eig.val))


    res.vec[i]=sum(eig.val)
  }


#  mej  <- cumsum(eig.val)
  Meff <- res.vec[length(res.vec)]

  return(list(Meff=Meff, mej=res.vec))
}

# Galwey method to estimate Meff & mej
Galwey.me <- function(sigma.matrix) {
  if (class(sigma.matrix) != "matrix") {
    stop("sigma.matrix argument should be a numeric matrix.")
  }

  N <- ncol(sigma.matrix)

  f <- function(x) { ifelse(x >= 1, 1, 0) + (x - floor(x)) }

  res.vec=rep(NA,times=N)
  res.vec[1]=1
  for (i in 2:N){
  	eig = eigen(sigma.matrix[1:i,1:i])$values
  	eig <- eig[which(eig >= 0)]
  	res.vec[i] <- sum(sqrt(eig))^2/sum(eig)
  }


#  mej  <- cumsum(eig.val)
  Meff <- res.vec[length(res.vec)]

  return(list(Meff=Meff, mej=res.vec))
}
