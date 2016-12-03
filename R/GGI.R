GGI <- function(Y, snpX, genes.length = NULL, genes.info = NULL,
            method = c("minP","PCA", "CCA", "KCCA","CLD","PLSPM","GBIGM",
                       "GATES", "tTS", "tProd"), ...){

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  # Arguments checks
  if (nlevels(as.factor(Y)) != 2) {
    stop("response variable should be binary. (2 modes).")
  } else if (class(snpX) != "SnpMatrix") {
    stop("snpX argument should be SnpMatrix object.")
  } else if (length(Y) != nrow(snpX)) {
    stop("Response variable should be conformant with genes matrix rows number.")
  } else if (is.null(genes.length) && is.null(genes.info)) {
    stop("Genes must be defined with either genes.length or genes.info arguments.")
  } else if (!is.null(genes.length) && !is.numeric(genes.length)) {
    stop("genes.length should be a numeric vector.")
  } else if (!is.null(genes.length) && any(genes.length <= 0)) {
    stop("Genes length can't be lesser or equal to 0.")
  } else if (!is.null(genes.length) && length(genes.length) < 2) {
    stop("Select at least two genes.")
  } else if (!is.null(genes.length) &&  sum(genes.length) > ncol(snpX)) {
    stop("genes.length argument non conformant with snps.matrix rows count.")
  } else if (!is.null(genes.length) && any(is.na(genes.length))) {
    stop("When provided, genes.length argument can't have missing values (NA).")
  } else if (!is.null(genes.info) && (!(is.data.frame(genes.info) | nrow(genes.info) > ncol(snpX)))) {
    stop("When provided, genes.info should be a data.frame with less rows than or as much as snpX columns.")
  } else if (!is.null(genes.info) && ncol(genes.info) != 4) {
    stop("genes.info argument should have four columns. See help file.")
  } else if (!is.null(genes.info) && !all(names(genes.info) %in% c("Genenames", "SNPnames", "Position", "Chromosome"))) {
    stop("genes.info argument should have its columns named: Genenames, SNPnames, Position, Chromosome")
  } else if (!is.null(genes.info) && is.character(genes.info$Genenames)) {
    stop("gene.info argument's Gene.name column should be of class character.")
  } else if (!is.null(genes.info) && nlevels(genes.info$Genenames) < 2) {
    stop("Select at least two genes.")
  } else if (!is.null(genes.info) && is.character(genes.info$SNPnames)) {
    stop("gene.info argument's SNP.name column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$Position)) {
    stop("gene.info argument's Position column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$Chromosome)) {
    stop("gene.info argument's Chromosome column should be of class character.")
  } else if (!is.null(genes.info) && any(is.na(genes.info))) {
    stop("When provided, genes.info can't have missing values (NA).")
  } else if (any(is.na(snpX))) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (any(is.na(Y))) {
    stop("The response variable must be complete. No NAs are allowed.")
  } else if (!is.null(genes.length) && !is.null(genes.info)) {
    warning("Both genes.length and genes.info arguments were provided, only genes.info will be used.")
  }

  method <- match.arg(method)

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  #SnpMatrix and dataframe are reordered to make sure that all SNP of a gene are contiguous.
  if (!is.null(genes.info)){
    genes.info <- genes.info[order(genes.info$Genenames, genes.info$SNPnames), ]
    snpX <- snpX[, as.character(genes.info$SNPnames)]
  }

  #Genes names computing
  if (!is.null(genes.info)) {
    genes.names <- levels(genes.info$Genenames)
  } else if (!is.null(genes.length) && is.null(names(genes.length))){
    genes.names <- paste("Gene", seq(1:length(genes.length)), sep=".")
    names(genes.length) <- genes.names
  } else {
    genes.names <- names(genes.length)
  }

  #Interactions listing
  interactions <- combn(genes.names, m=2)

  #Indexes of the genes
  if (!is.null(genes.info)) {
    gene.start <- NULL
    gene.end <- NULL
    for (i in 1:nlevels(genes.info$Genenames)) {
      gene <- genes.info[which(genes.info$Genenames %in% levels(genes.info$Genenames)[i]), ]
      gene.start <- c(gene.start, min(which(colnames(snpX) %in% gene$SNPnames), na.rm = TRUE))
      gene.end   <- c(gene.end, max(which(colnames(snpX) %in% gene$SNPnames), na.rm = TRUE))
    }
  } else {
    gene.start <- c(0, cumsum(genes.length)[1:(length(genes.length) - 1)]) + 1
    gene.end   <- cumsum(genes.length)
  }

  names(gene.start) <- genes.names
  names(gene.end)   <- genes.names

  #Setup of the return object
  pval.matrix <- diag(0, length(genes.names))
  colnames(pval.matrix) <- genes.names
  rownames(pval.matrix) <- genes.names

  stat.matrix <- diag(0, length(genes.names))
  colnames(stat.matrix) <- genes.names
  rownames(stat.matrix) <- genes.names

	if (method=="PCA"){
		df.matrix <- diag(0, length(genes.names))
		colnames(df.matrix) <- genes.names
		rownames(df.matrix) <- genes.names
	}

	res.method <- method
	res.parameter <- NULL

  #Application of the method on the interactions
  for (i in 1:ncol(interactions)) {
    print(paste("Interaction between", interactions[1, i], "&", interactions[2, i], "-",i,"/",ncol(interactions)))

    G1 <- snpX[, gene.start[interactions[1, i]]:gene.end[interactions[1, i]]]
    G2 <- snpX[, gene.start[interactions[2, i]]:gene.end[interactions[2, i]]]

#    if(!method %in% c("minP","GATES","tTS","tProd") || ncol(G1)*ncol(G2)<1000){
    	tmp <- switch(method,
                                                                         CCA = CCA.test(Y, G1, G2, ...),
                                                                         KCCA = KCCA.test(Y, G1, G2, ...),
                                                                         CLD = CLD.test(Y, G1, G2, ...),
                                                                         PLSPM = PLSPM.test(Y, G1, G2, ...),
                                                                         GBIGM = GBIGM.test(Y, G1, G2, ...),
                                                                         PCA = PCA.test(Y, G1, G2, ...),
                                                                         minP = minP.test(Y, G1, G2, ...),
                                                                         GATES = gates.test(Y, G1, G2, ...),
                                                                         tTS   = tTS.test(Y, G1, G2, ...),
                                                                         tProd = tProd.test(Y, G1, G2, ...))
    pval.matrix[interactions[1, i], interactions[2, i]] <- tmp$p.value
    stat.matrix[interactions[1, i], interactions[2, i]] <- tmp$statistic
    if (method=="PCA"){
    	df.matrix[interactions[1, i], interactions[2, i]] <- tmp$parameters["df"]    	
    }
#    if (is.null(res.method)){res.method <- tmp$method}
    if (is.null(res.parameter)){res.parameter <- tmp$parameter}

#} 
#else {
#      warning("Too much interactions to test for SSI method (>1000). NA returned")
#      pval.matrix[interactions[1, i], interactions[2, i]] <- NA
#      stat.matrix[interactions[1, i], interactions[2, i]] <- NA
#    }
  }

#  genes.interactions[lower.tri(genes.interactions)] <- t(genes.interactions)[lower.tri(genes.interactions)]
  pval.matrix[lower.tri(pval.matrix)] <- t(pval.matrix)[lower.tri(pval.matrix)]
  stat.matrix[lower.tri(stat.matrix)] <- t(stat.matrix)[lower.tri(stat.matrix)]
  if (method=="PCA"){
  	df.matrix[lower.tri(df.matrix)] <- t(df.matrix)[lower.tri(df.matrix)]
  	res <- list(statistic=stat.matrix,p.value=pval.matrix,df=df.matrix,method=res.method,parameter=res.parameter)
  	} else {
  		res <- list(statistic=stat.matrix,p.value=pval.matrix,method=res.method,parameter=res.parameter)
  	}
  	class(res) <- "GGInetwork"

  	return(res)
}

# Function that order genes.info and SnpMatrix so that all SNP of a genes are contiguous
order.snpMatrix <- function(snpX, genes.info) {
  genes.info <- genes.info[order(genes.info$Genenames, genes.info$SNPnames), ]
  snpX <- snpX[, as.character(genes.info$SNPnames)]


  return(list(snpX = snpX, genes.info = genes.info))
}
