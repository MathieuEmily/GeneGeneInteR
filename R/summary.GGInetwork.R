summary.GGInetwork <- function(object, ...){
	op <- options()
	options(scipen=-3,digits=2)
	cat("Gene-gene interaction network of ",ncol(object$p.value)," genes performed with:\n \t",object$method,"\n" )
	
	nn <- row.names(object$p.value)
	pval.none <- unlist(sapply(1:(length(nn)-1),FUN=function(i){object$p.value[i,(i+1):ncol(object$p.value)]}))

	pval <- data.frame(
	G1=unlist(sapply(1:(length(nn)-1),FUN=function(i){rep(nn[i],times=(length(nn)-i))})),
	G2=unlist(sapply(1:(length(nn)-1),FUN=function(i){nn[(i+1):length(nn)]})
),
	None=pval.none,
	bonferroni=p.adjust(pval.none,method="bonferroni"),
	BH=p.adjust(pval.none,method="BH")
	)
	
	w <- which(pval$None < 0.05)
	if (length(w)){
		tmp.df <- pval[w,c(1,2,3)]
		tmp.df <- tmp.df[order(tmp.df[,3]),]
		row.names(tmp.df) <- NULL
		names(tmp.df) <- c("Gene1","Gene2","Uncorrected p-value")
		cat("\nSignificant interaction with no correction at the level of 0.05 \n-------\n")
		print(tmp.df)
	} else{
		cat("\nNo significant interaction (at the level of 0.05) with no correction\n")
		}
	w <- which(pval$bonferroni < 0.05)
	if (length(w) > 0){
		tmp.df <- pval[w,c(1,2,4)]
		tmp.df <- tmp.df[order(tmp.df[,3]),]
		row.names(tmp.df) <- NULL
		names(tmp.df) <- c("Gene1","Gene2","bonferroni p-value")
		cat("\nSignificant interaction with a bonferroni correction at the level of 0.05 \n-------\n")
		print(tmp.df)
	} else {
		cat("\nNo significant interaction (at the level of 0.05) with a bonferroni correction\n")
	}
	w <- which(pval$BH < 0.05)
	if (length(w) > 0){
		tmp.df <- pval[w,c(1,2,5)]
		tmp.df <- tmp.df[order(tmp.df[,3]),]
		row.names(tmp.df) <- NULL
		names(tmp.df) <- c("Gene1","Gene2","BH p-value")
		cat("\nSignificant interaction with a Benjamini & Hochberg correction at the level of 0.05 \n-------\n")
		print(tmp.df)
	} else {
		cat("\nNo significant interaction (at the level of 0.05) with a Benjamini & Hochberg correction\n")
	}
	options(scipen=op$scipen,digits=op$digits)
}
