snpMatrixScour <- function(snpX, genes.length = NULL, genes.info = NULL,
                           min.maf = 0.01, min.eq = 0.01, call.rate=0.9) {
  if (class(snpX) != "SnpMatrix") {
    stop("snpX argument should be SnpMatrix object.")
  } else if (!is.null(genes.length) && (!(is.numeric(genes.length) || sum(genes.length) > ncol(snpX)))) {
    stop("genes.length argument should be a numeric vector which sum should be lesser than ncol(snpX).")
  } else if (!is.null(genes.length) && any(is.na(genes.length))) {
    stop("genes.length argument can't have missing values (NA).")
  } else if (!is.null(genes.info) && (!(is.data.frame(genes.info) || nrow(genes.info) > ncol(snpX)))) {
    stop("genes.info should be a data.frame with less rows than or as much as snpX columns.")
  } else if (!is.null(genes.info) && ncol(genes.info) != 4) {
    stop("genes.info argument should have four columns.")
  } else if (!is.null(genes.info) && !all(names(genes.info) %in% c("Genenames", "SNPnames", "Position", "Chromosome"))) {
    stop("genes.info argument should have its columns named: Genenames, SNPnames, Position, Chromosome")
  } else if (!is.null(genes.info) && is.character(genes.info$Genenames)) {
    stop("gene.info argument's Gene.name column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$SNPnames)) {
    stop("gene.info argument's SNP.name column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$Position)) {
    stop("gene.info argument's Position column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$Chromosome)) {
    stop("gene.info argument's Chromosome column should be of class character.")
  } else if (!is.null(genes.info) && any(is.na(genes.info))) {
    stop("genes.info can't have missing values (NA).")
  } else if (!is.numeric(min.maf) || min.maf < 0 || min.maf > 0.5) {
    stop("min.maf argument should be a numeric between 0 & 0.5.")
  } else if (!is.numeric(min.eq) || min.eq < 0 || min.eq > 1) {
    stop("min.eq argment should be a numeric between 0 & 1.")
  } else if (!is.numeric(call.rate) || call.rate > 1 || call.rate < 0) {
    stop("call.rate argment should be a numeric between 0 & 1.")
  } else if (!is.null(genes.length) && !is.null(genes.info)) {
    stop("Both genes.length and genes.info arguments were provided, only genes.info will be used.")
  }

  # If no genes were defined, whole dataset is scoured as is
  if (is.null(genes.length) && is.null(genes.info)) {
  	
    new.snpX <- GeneScour(snpX, min.maf, min.eq, call.rate)
  	new.snpX <- as(new.snpX, "SnpMatrix")
	 return(list(snpX = new.snpX, genes.info = NULL))

  # If genes are defined, each genes is scoured individually to
  # keep track of its length variation.
  } else {
    #Indexes of the genes
    if (!is.null(genes.info)) {

      #SnpMatrix and dataframe are reordered to make sure that all SNP of a gene are contiguous.
      genes.info <- genes.info[order(genes.info$Genenames, genes.info$SNPnames), ]
      snpX <- snpX[, as.character(genes.info$SNPnames)]

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

    n.genes <- length(gene.start)
    new.snpX <- NULL
    for (i in 1:n.genes) {
      cur.gene <- try(GeneScour(snpX[, gene.start[i]:gene.end[i]], min.maf, min.eq, call.rate), silent=TRUE)

      if (class(cur.gene) == "try-error") {
        warning("A whole gene had to be removed as no SNP met the requirements.")
      }

      # Gene positions are updated
      if (is.null(genes.info)) {
        if (class(cur.gene) == "try-error"){
          new.genes[i] <- NA
        } else {
          new.genes[i] <- ncol(cur.gene)
        }
      } else {
        if (class(cur.gene) == "try-error"){
          cur.gene.name <- levels(genes.info$Genenames)[i]
          genes.info <- genes.info[- which(genes.info$Genenames == cur.gene.name),]
        } else {
          SNP.cond  <- !(genes.info$SNPnames %in% colnames(cur.gene))
          gene.cond <- genes.info$Genenames == levels(genes.info$Genenames)[i]
          select <- which(SNP.cond & gene.cond)
          if (length(select) > 0) {genes.info <- genes.info[- select,]}
        }
      }

      # Binding matrices
      if (class(cur.gene) != "try-error" && is.null(new.snpX)) {
        new.snpX <- cur.gene
      } else if (class(cur.gene) != "try-error") {
        new.snpX <- cbind(new.snpX, cur.gene)
      }
    }

    if (!is.null(genes.info)) {
      new.genes <- droplevels(genes.info)
    } else if (length(which(is.na(new.genes))) > 0) {
      new.genes <- new.genes[-which(is.na(new.genes))]
    }
  }

  new.snpX <- as(new.snpX, "SnpMatrix")

  print("A list object has been returned with elements: snpX & genes.info")
  return(list(snpX = new.snpX, genes.info = new.genes))
}

# Function that checks SNP validity
# A SNP is considered valid if it meets following criteria:
#     - MAF is superior to min.maf arguments
#     - Hardy-Weinberg equilibrium is verified using min.eq threshold
GeneScour <- function(gene, min.maf = 0.01, min.eq = 0.01, call.rate=0.9){
  # NA filter
#  SNP.CallRate <- colSums(!is.na(gene))/nrow(gene)
	SNP.CallRate <- snpStats::col.summary(gene)$Call.rate
	
  # MAF filtering
  SNP.MAF <- snpStats::col.summary(gene)$MAF

  # Hardy-Weinberg Equilibrium check
  # AA <- colSums(as(gene, "numeric") == 2, na.rm=TRUE)
  # Aa <- colSums(as(gene, "numeric") == 1, na.rm=TRUE)
  # aa <- colSums(as(gene, "numeric") == 0, na.rm=TRUE)

  # p <- (2*AA + Aa)/(2*(AA + Aa + aa))
  # q <- 1 - p

  # E.AA <- p^2 * nrow(gene)
  # E.Aa <- 2 * p * q * nrow(gene)
  # E.aa <- q^2 * nrow(gene)

  # # Chi-Square test
  # X.AA <- ((AA - E.AA)^2)/E.AA
  # X.Aa <- ((Aa - E.Aa)^2)/E.Aa
  # X.aa <- ((aa - E.aa)^2)/E.aa
  # X <- X.AA + X.Aa + X.aa
  # p.X <- pchisq(X, 1, lower.tail=FALSE)
	
	p.X <- 2*(1-pnorm(abs(snpStats::col.summary(gene)$z.HWE)))


  # Filtering
  if (any(SNP.MAF >= min.maf & p.X >= min.eq & SNP.CallRate >= call.rate)) {
    return(gene[, (SNP.MAF >= min.maf & p.X >= min.eq & SNP.CallRate >= call.rate)])
  } else {
    stop("No SNP with sufficient MAF and verifying Hardy-Weinberg equilibrium")
  }
}
