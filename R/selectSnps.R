selectSnps <- function(snpX, genes.info, select){
  if(class(snpX)!="SnpMatrix"){
    stop("snpX must be a snpMatrix object.")
  } else if(class(genes.info)!="data.frame"){
    stop("genes.info must be a data.frame.")
  } else if(ncol(snpX)!=nrow(genes.info)){
    stop("snpX and genes.info must have the same number of snps.")
  }

  res <- NULL

  if(class(select)%in%c("numeric","integer")){
    if(any(c(select<1, select>ncol(snpX)))){
      stop("The selection of snps is out of bounds.")
    } else {
      res[["snpX"]] <- snpX[,select]
      res[["genes.info"]] <- genes.info[select,]
    }

  } else if(class(select)=="character"){
    if(all(select %in% colnames(snpX))){
      res[["snpX"]] <- snpX[,select]
      res[["genes.info"]] <- genes.info[genes.info[,"SNPnames"]%in%select,]
    } else if(all(select %in% genes.info[,"Genenames"])){
      genes <- genes.info[,"Genenames"]%in%select
      res[["snpX"]] <- snpX[,as.character(genes.info[genes,"SNPnames"])]
      res[["genes.info"]] <- genes.info[genes,]
    } else {
      liste <- strsplit(select, ":")
      snp <- c()
      if(length(liste[[1]])==2){
        for(i in 1:length(liste)){
          rows <- genes.info[as.numeric(genes.info[,"Position"])>=as.numeric(liste[[i]][1]),]
          rows <- rows[as.numeric(rows[,"Position"])<=as.numeric(liste[[i]][2]),]
          snp <- c(snp,as.character(rows[,"SNPnames"]))
        }
      } else if(length(liste[[1]])==3){
        for(i in 1:length(liste)){
          rows <- genes.info[as.numeric(genes.info[,"Chromosome"])==as.numeric(liste[[i]][1]),]
          rows <- rows[as.numeric(rows[,"Position"])>=as.numeric(liste[[i]][2]),]
          rows <- rows[as.numeric(rows[,"Position"])<=as.numeric(liste[[i]][3]),]
          snp <- c(snp,as.character(rows[,"SNPnames"]))
        }
      } else {stop("Wrong format for select argument, or genes or snps not found. Please refer to help file.")}

      if(length(snp)==0){stop("No matches were found with your selection.")}
      else{
        res[["snpX"]] <- snpX[,snp]
        res[["genes.info"]] <- genes.info[genes.info[,"SNPnames"]%in%snp,]
      }
    }

  } else {
    stop("The argument select must be a numeric or character vector.")
  }

  rownames(res[["genes.info"]]) <- 1:nrow(res[["genes.info"]])
  res$genes.info <- droplevels(res$genes.info)
  return(res)
}
