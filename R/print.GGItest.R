print.GGItest <- function(x, ...){
	cat("Gene-Gene Interaction method performed with:\n \t",x$method,"\n")
	if (is.null(x$df)){
		cat(names(x$statistic)," = ",x$statistic,", p-value = ",x$p.value,"\n",sep="")
	} else {
		cat(names(x$statistic)," = ",x$statistic,", df = ",x$df,", p-value = ",x$p.value,"\n",sep="")	
		}
}




