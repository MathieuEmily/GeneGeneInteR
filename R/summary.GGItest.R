summary.GGItest <- function(object, ...){
    x <- object
	cat("Gene-Gene Interaction method performed with:\n \t",x$method,"\n")
	if (is.null(x$df)){
		cat(names(x$statistic)," = ",x$statistic,", p-value = ",x$p.value,"\n",sep="")
	} else {
		cat(names(x$statistic)," = ",x$statistic,", df = ",x$df,", p-value = ",x$p.value,"\n",sep="")	
		}
	ll <- x$parameter
	if (!is.null(x$parameter)){
		cat("-- ")
		cat("List of parameter(s):\n--- ")
		for (i in 1:length(ll)){
			 cat(names(ll)[i], "=", ll[[i]])
            if (i != length(ll)) cat(", ")
            if (i == length(ll)) cat("\n ")
		}
	}
}

