PCA.test <- function(Y,G1,G2,threshold=0.8,method="GenFreq"){
	if (method=="GenFreq"){
		tmp <- PCA.GenFreq(Y=Y,G1=G1,G2=G2,threshold=threshold)
		tmp$data.name <- paste(deparse(substitute(Y))," and  (",deparse(substitute(G1))," , ",deparse(substitute(G2)),")",sep="")
		tmp$method <- paste(tmp$method,"-",method)
		tmp$parameters <- c(tmp$parameters,threshold)
		names(tmp$parameters) <- c("df","threshold")
		return(tmp)
	} else if (method=="Std"){
		tmp <- PCA.Std(Y=Y,G1=G1,G2=G2,threshold=threshold)
		tmp$data.name <- paste(deparse(substitute(Y))," and  (",deparse(substitute(G1))," , ",deparse(substitute(G2)),")",sep="")
		tmp$method <- paste(tmp$method,"-",method)
		tmp$parameters <- c(tmp$parameters,threshold)
		names(tmp$parameters) <- c("df","threshold")
		return(tmp)
	} else {
		stop("method argument should be a character string either GenFreq or Std")
	}
}


PCA.Std <- function(Y, G1, G2, threshold=0.8) {

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  # Arguments checks
  if (class(threshold) != "numeric") {
    stop("threshold argument should be a numeric.")
  } else if (threshold < 0 | threshold > 1) {
    stop("threshold argument shoud be comprised in [0, 1] interval.")
  } else if (nlevels(as.factor(Y)) != 2) {
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

  # Threshold formatting
  inertia.thresh <- threshold * 100

  # SnpMatrix coerced into matrix to be compatible with FactoMineR::PCA
  G1.num <- as(G1, "numeric")
  G2.num <- as(G2, "numeric")

  # PCA performed
  G1.PCA <- FactoMineR::PCA(G1.num, ncp=NULL, graph=FALSE)
  G2.PCA <- FactoMineR::PCA(G2.num, ncp=NULL, graph=FALSE)

  # Genes are represented by PCA coords on enough dimensions to retrieve
  # as much inertia as set by user.
  G1.ncp <- which(G1.PCA$eig[, 3] > inertia.thresh)[1]
  G2.ncp <- which(G2.PCA$eig[, 3] > inertia.thresh)[1]

  if (G1.ncp == 1) {
    G1.VarSynth <- data.frame(Dim.1 =G1.PCA$ind$coord[, 1:G1.ncp])
  } else {
    G1.VarSynth <- G1.PCA$ind$coord[, 1:G1.ncp]
  }

  if (G2.ncp == 1) {
    G2.VarSynth <- data.frame(Dim.1 =G2.PCA$ind$coord[, 1:G2.ncp])
  } else {
    G2.VarSynth <- G2.PCA$ind$coord[, 1:G2.ncp]
  }

  G1.PCA <- list(VarSynth = G1.VarSynth, Inertia = G1.PCA$eig[1:G1.ncp, 2])
  G2.PCA <- list(VarSynth = G2.VarSynth, Inertia = G2.PCA$eig[1:G2.ncp, 2])

  # Interaction effects are tested
  return(compare.PCA(Y, G1.PCA, G2.PCA))
}

#'@describeIn PCA.Std Standardization based on Hardy-Weinberg equilirum
#'@export
PCA.GenFreq <- function(Y, G1, G2, threshold=0.8) {

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  # Arguments checks
  if (class(threshold) != "numeric") {
    stop("thresold argument should be a numeric.")
  } else if (threshold < 0 | threshold > 1) {
    stop("threshold argument shoud be comprised in [0, 1] interval.")
  } else if (nlevels(as.factor(Y)) != 2) {
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

  # Threshold formatting
  inertia.thresh <- threshold * 100

  G1.PCA <- get.PCA.res(G1)
  G2.PCA <- get.PCA.res(G2)

  # Genes are represented by PCA coords on enough dimensions to retrieve
  # as much inertia as set by user.
  G1.ncp <- which(cumsum(G1.PCA$Inertia) > inertia.thresh)[1]
  G2.ncp <- which(cumsum(G2.PCA$Inertia) > inertia.thresh)[1]

  G1.PCA <- list(VarSynth=G1.PCA$VarSynth[, 1:G1.ncp], Inertia=G1.PCA$Inertia[1:G1.ncp])
  G2.PCA <- list(VarSynth=G2.PCA$VarSynth[, 1:G2.ncp], Inertia=G2.PCA$Inertia[1:G2.ncp])

  if (G1.ncp == 1) {
    G1.PCA$VarSynth <- data.frame(Dim.1 = G1.PCA$VarSynth)
  }

  if (G2.ncp == 1) {
    G2.PCA$VarSynth <- data.frame(Dim.1 = G2.PCA$VarSynth)
  }

  # Interaction effects are tested
  return(compare.PCA(Y, G1.PCA, G2.PCA))
}

## Function that retrieves Principal Components and Eigen Values
## of a SnpMatrix object.
get.PCA.res <- function(gene.matrix){
  if (class(gene.matrix) != "SnpMatrix") {
    stop("gene.matrix argument should be SnpMatrix object.")
  } else if (sum(is.na(gene.matrix))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  }

  Id <- diag(ncol(gene.matrix))
  # X'X matrix is calculated
  G.XpX <- snpStats::snp.pre.multiply(gene.matrix, t(snpStats::snp.post.multiply(gene.matrix, Id)))
  # Singular Values Decomposition is performed
  G.eigen <- eigen(G.XpX, symmetric=TRUE)
  # Individuals coordinates are calculated
  G.PCA <- snpStats::snp.post.multiply(gene.matrix, G.eigen$vectors)
  # Columns names are formatted (for comparison procedure)
  colnames(G.PCA) <- paste("Dim", seq(1, ncol(G.PCA)), sep=".")

  # Computes inertia percentages for all ordered eigen values
  eigen.val <- 100*G.eigen$values/sum(G.eigen$values)

  return(list(VarSynth=G.PCA, Inertia=eigen.val))
}

## Function that fits single effects model and interaction
## model and compare them.
compare.PCA <- function(Resp., G1.PCA, G2.PCA) {
  # Arguments checks
  if (nlevels(as.factor(Resp.)) != 2) {
    stop("response variable should be binary. (2 modes).")
  } else if (class(G1.PCA) != "list" | class(G2.PCA) != "list") {
    stop("G1.PCA and G2.PCA arguments should be list objects.")
  } else if (length(G1.PCA) != 2 | length(G2.PCA) != 2) {
    stop("G1.PCA and G2.PCA arguments should be of length 2: VarSynth and Inertia elements.")
  } else if (nrow(G1.PCA$VarSynth) != nrow(G2.PCA$VarSynth)) {
    stop("VarSynth elements of G1.PCA and G2.PCA should have same rows count.")
  } else if (length(Resp.) != nrow(G1.PCA$VarSynth) | length(Resp.) != nrow(G2.PCA$VarSynth)) {
    stop("Response variable should be conformant with genes matrices rows number.")
  } else if (length(G1.PCA$Inertia) > ncol(G1.PCA$VarSynth) | length(G2.PCA$Inertia) > ncol(G2.PCA$VarSynth)){
    stop("Inertia elements of G1.PCA and G2.PCA can't be of greater length than ncol(VarSynth)")
  } else if (sum(is.na(Resp.)) != 0) {
    stop("The response variable vector must be complete. No NAs are allowed.")
  }

  G1.inertia <- G1.PCA$Inertia
  G2.inertia <- G2.PCA$Inertia
  G1.PCA <- G1.PCA$VarSynth
  G2.PCA <- G2.PCA$VarSynth

  # Effects' names
  colnames(G1.PCA) <- paste("G1", colnames(G1.PCA), sep=".")
  colnames(G2.PCA) <- paste("G2", colnames(G2.PCA), sep=".")

  # Iterating until model is fitted
  # Less informant component is removed each time fitting is impossible.
  keep.on <- TRUE
  trimmed <- FALSE
  while (keep.on) {
    data <- data.frame(Resp., G1.PCA, G2.PCA)

    # Effects set up
    single.effects <- paste(colnames(G1.PCA), colnames(G2.PCA), sep="+", collapse="+")

    # Creating all possible pairs
    inter.effects <- expand.grid(colnames(G1.PCA), colnames(G2.PCA))
    # Creating a string with all interactions
    inter.effects <- paste(inter.effects[, 1], inter.effects[, 2], sep=":", collapse="+")

    inter.formula <- as.formula(paste("Resp.~", single.effects, "+", inter.effects, sep=""))
    null.formula <- as.formula(paste("Resp.", single.effects, sep="~"))

    # Trying to fit the largest model
    inter.mod <- tryCatch(glm(inter.formula, data=data, family = "binomial"),
             warning = function(w){"warning"},
             error   = function(e){"error"})

    # If succeeded end the process
    if (any(class(inter.mod) == "glm")){
      keep.on <- FALSE
      null.mod  <- glm(null.formula, data=data, family = "binomial")
      comp.res <- anova(null.mod, inter.mod, test="Chisq")
    # If model coulnd't be fitted then a component is removed
    } else {
      trimmed <- TRUE
      # Looking for the less informant component
      if (G1.inertia[length(G1.inertia)] < G2.inertia[length(G2.inertia)]) {
		# If less informant component is one of the two first PC of the gene
		# then a PC is removed from the other gene instead.
        if (length(G1.inertia) > 2) {
          G1.inertia <- G1.inertia[-length(G1.inertia)]
          G1.PCA <- G1.PCA[, -ncol(G1.PCA)]
        } else if (length(G2.inertia) > 2){
          G2.inertia <- G2.inertia[-length(G2.inertia)]
          G2.PCA <- G2.PCA[, -ncol(G2.PCA)]
		# If both genes are only described by their first two PC then an error is issued.
        } else {
          stop("Genes too correlated to fit glm model.")
        }
      } else {
        if (length(G2.inertia) > 2) {
          G2.inertia <- G2.inertia[-length(G2.inertia)]
          G2.PCA <- G2.PCA[, -ncol(G2.PCA)]
        } else if (length(G1.inertia) > 2){
          G1.inertia <- G1.inertia[-length(G1.inertia)]
          G1.PCA <- G1.PCA[, -ncol(G1.PCA)]
        } else {
          stop("Genes too correlated to fit glm model.")
        }
      }
    }
  }

  # When PC had to be removed, a warning is issued to inform the user of
  # inertia loss.
  if (trimmed) {
    warn.str <- paste("Less principal components had to be kept to fit glm model: ",
                      round(max(cumsum(G1.inertia)),2), "% of inertia was kept for G1 & ",
                      round(max(cumsum(G2.inertia)),2), "% of inertia was kept for G2.", sep="")
    warning(warn.str)
  }

pval <- comp.res$'Pr(>Chi)'[2]
stat <- comp.res$'Deviance'[2]
df <- comp.res$'Df'[2]
names(stat)="Deviance"
null.value <- 0
names(null.value) <- "deviance"

estimate <- c(comp.res$"Resid. Dev"[1],comp.res$"Resid. Dev"[2])
names(estimate) <- c("Deviance without interaction","Deviance with interaction")
	res <- list(
		null.value=null.value,
		alternative="greater",
		method="Gene-based interaction based on Principal Component Analysis",
		estimate= estimate,
		data.name="d",
		statistic=stat,
		p.value=pval,
		parameters=df)
	class(res) <- "htest"
  return(res)
#  return(comp.res$'Pr(>Chi)'[2])
}
