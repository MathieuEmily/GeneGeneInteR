chclust <- function(d, method="coniss") {
  if (!("dist" %in% class(d)))
     stop("Input must be a distance matrix")
  x <- as.matrix(d)
  if (!is.numeric(d))
     stop("Input matrix must be numeric")
  if (any(is.na(d)))
     stop("Missing values in input data")
  METHODS <- c("conslink", "coniss")
  method <- pmatch(method, METHODS)
  if(is.na(method))
     stop("Invalid clustering method")
  if(method == -1)
     stop("Ambiguous clustering method")
 ret <- .Call("chclust_c", x, as.integer(method), PACKAGE="GeneGeneInteR")
 if (is.character(ret))
     stop(ret)
  merge <- .find.groups(ret)
  tree <- list(merge=merge[, seq_len(2)], height=sort(ret), seqdist = merge[, 3], order=seq_len(nrow(x)), labels=attr(d, "Labels"), method=METHODS[method], call=match.call(), dist.method = attr(d, "method"))
  class(tree) <- c("chclust", "hclust")
  tree
}

.find.groups <- function(height) {
  nr = length(height)
  x <- height
  merge <- matrix(nrow=nr, ncol=3)
  rec <- vector(mode="numeric", length=nr+1) 
  rec[] <- NA
  nG = 1
  for (i in seq_len(nr)) {
     n <- which.min(x)
     minx <- min(x, na.rm=TRUE)
     merge[i, 3] <- minx
     if (is.na(rec[n]) & is.na(rec[n+1])) {
        merge[i,1] = -n
        merge[i,2] = -(n+1)
        rec[n] = nG
        rec[n+1] = nG
     } else {
        if (is.na(rec[n]) & !is.na(rec[n+1])) {
           merge[i,1] = -n
           merge[i,2] = rec[n+1]
           rec[n] = nG
           rec[rec == rec[n+1]] = nG
        } else {
           if (!is.na(rec[n]) & is.na(rec[n+1])) {
              merge[i,1] = rec[n]
              merge[i,2] = -(n+1)
              rec[rec == rec[n]] = nG
              rec[n+1] = nG
           } else {
              merge[i,1] = rec[n]
              merge[i,2] = rec[n+1]
              rec[rec == rec[n]] = nG
              rec[rec == rec[n+1]] = nG
           }
        }
     }
     x[n] <- NA
     nG <- nG+1
  }
  merge
}

