importFile <- function (file, pos, pos.sep="\t", ...) {

  if(class(file)!="character"){
    stop("Please, enter a valid path file for genotype data.")
  }

  nc <- nchar(file)
  ext <- substr(file, nc - 3, nc)
  extgz <- substr(file, nc - 6, nc)
  res <- list()
  imp <- NULL

  if (ext==".ped") {
    imp <- snpStats::read.pedfile(file = file, ...)
	res[["status"]]	<- as.factor(imp$fam$affected)

  } else if (extgz==".ped.gz") {
    imp <- snpStats::read.pedfile(file = file, ...)
	res[["status"]]	<- as.factor(imp$fam$affected)

  } else if (ext == ".bed") {
    imp <- snpStats::read.plink(bed = file, ...)

  } else if (ext %in% c(".vcf", "f.gz")) {
    imp <- GGtools::vcf2sm(tbxfi = Rsamtools::TabixFile(file), ...)

  } else if (ext == ".impute2") {
    imp <- snpStats::read.impute(file, ...)

  } else {stop("Please enter a valid pedfile, plink, vcf, or impute2 file.")}

  if(class(imp)=="list"){
    res[["snpX"]] <- imp$genotypes

  } else if(class(imp)=="SnpMatrix"){
    res[["snpX"]] <- imp
  }

  if(!missing(pos)){

    if(length(pos)==1){
      if(class(pos)=="character"){
        infos <- read.csv(pos, sep=pos.sep, header=TRUE)

        chr <- infos[,names(infos)%in%c("Chromosome","Chr","chromosome","chr")]
        gene <- infos[,names(infos)%in%c("Gene","gene","genenames","Genenames","Gene.names","gene.names")]
        snp <- infos[,names(infos)%in%c("SNP","Snp","snp","SNPnames","Snpnames","snpnames","SNP.names","Snp.names","snp.names")]
        posi <- infos[,names(infos)%in%c("Position","position","pos","Pos")]

        if(length(chr)==0){chr<-rep(NA,nrow(infos));warning("Chromosome column was not found.")}
        else if(length(gene)==0){gene<-rep(NA,nrow(infos));warning("Gene names column was not found.")}
        else if(length(snp)==0){snp<-colnames(res[["snpX"]])}
        else if(length(posi)==0){posi<-rep(NA,nrow(infos));warning("Position column was not found.")}

        res[["genes.info"]] <- data.frame(Chromosome=chr,
                                    Genenames=gene,
                                    SNPnames=snp,
                                    Position=posi)
      } else {
        stop("Pos argument needs to be either a numeric vector, a character vector or a path file.")
      }
    } else if(length(pos)==ncol(res[["snpX"]])){
      if(class(pos)%in%c("numeric","integer")){
        snp <- names(pos)
        if(is.null(snp)){snp <- colnames(res[["snpX"]])}
        res[["genes.info"]] <- data.frame(Chromosome=rep(NA,length(pos)),
                                    Genenames=rep(NA,length(pos)),
                                    SNPnames=snp,
                                    Position=pos)

      } else if(class(pos)=="character"){
        snp <- names(pos)
        if(is.null(snp)){snp <- colnames(res[["snpX"]])}

        liste <- data.table::tstrsplit(pos, split=":", fixed=TRUE)
        res[["genes.info"]] <- data.frame(Chromosome=liste[[1]],
                                    Genenames=rep(NA,length(pos)),
                                    SNPnames=snp,
                                    Position=liste[[2]])
      } else {
        stop("Pos argument needs to be either a numeric vector, a character vector or a path file.")

      }
    } else {
      stop("The number of SNPs must be the same in genotype data and position information.")
    }
  } else {
    if(class(imp)=="list"){
      chr <- imp$map[,names(imp$map)%in%c("Chromosome","Chr","chromosome","chr")]
      gene <- imp$map[,names(imp$map)%in%c("Gene","gene","genenames","Genenames","Gene.names","gene.names")]
      snp <- imp$map[,names(imp$map)%in%c("SNP","Snp","snp","SNPnames","Snpnames","snpnames","SNP.names","Snp.names","snp.names","snp.name","SNP.name","Snp.name","SNPname","Snpname","snpname")]
      posi <- imp$map[,names(imp$map)%in%c("Position","position","pos","Pos","V2")]

      if(class(posi)=="character"){
        liste <- data.table::tstrsplit(posi, split=":", fixed=TRUE)
        chr <- liste[[1]]
        posi <- liste[[2]]
      }

      if(length(chr)==0){chr<-rep(NA,nrow(imp$map));warning("Chromosome column was not found.")}
      if(length(gene)==0){gene<-rep(NA,nrow(imp$map));warning("Gene names column was not found.")}
      if(length(snp)==0){snp<-colnames(res[["snpX"]])}
      if(length(posi)==0){posi<-rep(NA,nrow(imp$map));warning("Position column was not found.")}

      res[["genes.info"]] <- data.frame(Chromosome=chr,
                                  Genenames=gene,
                                  SNPnames=snp,
                                  Position=posi)
    } else {
      chr<-rep(NA,ncol(res[["snpX"]]))
      gene<-rep(NA,ncol(res[["snpX"]]))
      snp<-colnames(res[["snpX"]])
      posi<-rep(NA,ncol(res[["snpX"]]))

      res[["genes.info"]] <- data.frame(Chromosome=chr,
                                  Genenames=gene,
                                  SNPnames=snp,
                                  Position=posi)
    }
  }

  if(any(colnames(res[["snpX"]])!=res[["genes.info"]][,"SNPnames"])){warning("Be careful, the SNP names don't match between snpMatrix and info dataframe.")}

  return(res)
}
