imputeSnpMatrix <- function(snpX, genes.info,
                            on.rem = c("SNP", "ind", "none"), quiet=FALSE){
  if (class(snpX) != "SnpMatrix") {
    stop("snpX argument should be SnpMatrix object.")
  } else if (!is.null(genes.info) && (!(is.data.frame(genes.info) | nrow(genes.info) > ncol(snpX)))) {
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
  } else if (!is.logical(quiet)) {
    stop("quiet argument should be a logical.")
  }

  on.rem <- match.arg(on.rem)

  imputed <- as(snpX, "numeric")
  if (!quiet){prog <- txtProgressBar(0, nrow(snpX), char="-", style = 3)}
  for (i in 1:nrow(snpX)) {

    select <- which(is.na(snpX[i, ]))

    if (length(select) > 0) {
      missing <- snpX[-i, select]
      present <- snpX[-i, -select]

      pos.miss <- genes.info$Position[select]
      pos.pres <- genes.info$Position[-select]

      rules <- silence(snpStats::snp.imputation)(present, missing, pos.pres, pos.miss)

      imp.targ <- snpStats::impute.snps(rules, snpX[i, ])

      imputed[i, is.na(imputed[i, ])] <- round(imp.targ)
    }

    if (!quiet){setTxtProgressBar(prog, i)}
  }
  cat("\n")

  if (any(is.na(imputed)) && on.rem == "SNP") {
    select <- which(colSums(is.na(imputed)) > 0)
    genes.info <- droplevels(genes.info[!(genes.info$SNPnames %in% colnames(imputed)[select]), ])
    imputed <- imputed[, -select]

    warning(paste(length(select), "SNP were removed due to remaining missing values."))
  } else if (any(is.na(imputed)) && on.rem == "ind") {
    select <- which(colSums(is.na(imputed)) > 0)
    imputed <- imputed[-select,]

    warning(paste(length(select), "individuals were removed due to remaining missing values."))
  } else if (any(is.na(imputed))) {
    warning(paste(sum(is.na(imputed)), "remaining missing values."))
  }

  imputed <- as(imputed, "SnpMatrix")

  return(list(snpX = imputed, genes.info = genes.info))
}

# Function used to silence snpStats::snp.imputation
silence <- function(f){
  return(function(...) {capture.output(w<-f(...));return(w);});
}
