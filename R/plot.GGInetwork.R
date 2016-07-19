plot.GGInetwork <- function(x,method=c("heatmap","network"), threshold=NULL, col=c("#D6604D", "#104E8B"), colbar.width=0.15, title=NULL,  hclust.order=FALSE, use.log=FALSE, NA.col="#D3D3D3", draw.pvals=(ncol(x$p.value) <= 15), draw.names=(ncol(x$p.value) <= 25), interact=FALSE, method.adjust=c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr"), genes=1:ncol(x$p.value), plot.nointer=TRUE, ...){
	if(class(x)!="GGInetwork") {
	    stop("x should be an object of class GGInetwork.")
	}
	method <- match.arg(method)
	switch(method,
			heatmap=GGI.plot(GGI=x$p.value, genes=genes, col=col, colbar.width=colbar.width, title=title, hclust.order=hclust.order, use.log=use.log, threshold=threshold, NA.col=NA.col, draw.pvals=draw.pvals, interact=interact, method.adjust=method.adjust),
			network=draw.network(GGI=x$p.value, genes=genes, threshold=threshold, plot.nointer=plot.nointer))
}


GGI.plot <- function(GGI,genes=1:ncol(GGI),col=c("#D6604D", "#104E8B"), colbar.width=0.15,
                     title=NULL, hclust.order=FALSE, use.log=FALSE,
                     threshold=NULL, NA.col="#D3D3D3",
                     draw.pvals=(ncol(GGI) <= 15), draw.names=(ncol(GGI) <= 25),
                     interact=!(draw.pvals && draw.names),method.adjust=c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr")) {

  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.character(col)) {
    stop("col argument should be a character vector.")
  } else if (!is.numeric(colbar.width) || colbar.width < 0) {
    stop("colbar.width argument should be a positive numeric")
  } else if (!is.character(title) && !is.null(title)) {
    stop("title argument should be a string.")
  } else if (!is.logical(draw.pvals) | !is.logical(draw.names)) {
    stop("show.pvals & draw.names arguments should be logical.")
  } else if (!is.logical(hclust.order)) {
    stop("hclust.order argument should be logical.")
  } else if (!is.logical(use.log)) {
    stop("use.log argument should be logical")
  } else if (!is.logical(interact)) {
    stop("interact argument should be logical")
  } else if (!is.null(threshold) && is.numeric(threshold) && (threshold > 1 || threshold < 0)) {
    stop("threshold argument can not be a numeric greater than 1 or lesser than 0.")
  } else if (!is.null(threshold) && is.character(threshold) && threshold != "R") {
    stop("threshold argument can not be any other string than 'R'.")
  } else if (!is.character(NA.col)) {
    stop("NA.col argument should be a character.")
  }

  if(class(genes)=="character"&&any(!genes%in%colnames(GGI))){
    stop("Genes and GGI don't match. Please select genes that are named in GGI.")
  }

  GGI <- GGI[genes,genes]

	method.adjust <- match.arg(method.adjust)
	
	GGI[lower.tri(GGI)] <- p.adjust(GGI[lower.tri(GGI)],method=method.adjust)
	GGI[upper.tri(GGI)] <- p.adjust(GGI[upper.tri(GGI)],method=method.adjust)

  R.thresh <- c(0.001, 0.01, 0.05, 0.1)

  # If only one color is parsed, white is
  # used to complete the scale
  if (length(col) < 2) {
    col <- c(col, "#FFFFFF")
  }

  if (use.log){
    GGI <- -log10(GGI)
    diag(GGI) <- min(GGI[row(GGI) != col(GGI)])
    col <- rev(col)

    if (!is.null(threshold)) threshold <- -log10(threshold);

    R.thresh <- -log10(R.thresh)
  }

  col.FUN <- grDevices::colorRampPalette(col)

  # Names checking (generated if none)
  if (is.null(dimnames(GGI))){
    genes.names <- paste("Gene", 1:ncol(GGI), sep=".")
    dimnames(GGI) <- list(genes.names, genes.names)
  }

  # Clustering
  if (hclust.order) {
    GGI.clust <- hclust(as.dist(GGI))
    GGI <- GGI[GGI.clust$order, GGI.clust$order]
  }

  # Calculating plot size
  plot.setup(GGI, colbar.width, draw.names, threshold)

  # Draw color map
  rect.pos <- draw.matrix(GGI, col.FUN, threshold, NA.col, R.thresh, use.log)

  # Draw color legend bar
  leg <- draw.colbar(GGI, col.FUN, colbar.width, threshold, R.thresh, NA.col, use.log)

  if (!is.null(leg)) leg <- leg$rect

  # Draw genes names
  if (draw.names && ncol(GGI) <= 25) {
    draw.genes.names(dimnames(GGI), rect.pos)
  } else if (draw.names) {
    warning("GGI object is too big (26+ genes): genes names were not plotted.
            The use of the tooltip functionality is recommanded.")
  }

  # Draw p-values
  if (draw.pvals && ncol(GGI) <= 15) {
    draw.interp(GGI, rect.pos)
  } else if (draw.pvals) {
    warning("GGI object is too big (16+ genes): p-values were not plotted.
             The use of the tooltip functionality is recommanded.")
  }

  # Draw title
  if (is.null(title)) {
    title <- "Genes Interactions Matrix Plot"
  }
  title(main=title)

  # Activate tooltip functionality
  if (interact && interactive()) {
    writeLines("Click on a cell to get info on that cell.\nPress Esc. to leave.")

	  prime.plot <- recordPlot()
  	inter.tooltip  <- FALSE
  	keep.on    <- TRUE
  	while(keep.on) {
  	  coords <- locator(n=1)
  	  # If user forcefully stop the locator function then
  	  # break out of the loop.
  	  if (is.null(coords)) break
  	  # As the bottom left point is used to identify a square
  	  # on the plot, coordinates are floored.
  	  coords <- floor(as.numeric(coords))

  	  # Plot coordinates are converted back to matrix coordinates.
  	  coords <- c(row=nrow(GGI) - coords[2], col=coords[1])

  	  # Check if coordinates are conformant with GGI matrix
  	  coords.check <- try(GGI[coords['row'], coords['col']], silent=TRUE)
  	  if (class(coords.check) != 'try-error' && length(coords.check) == 1) {
  	    # Check if selected point is in upper triangle -diag excluded-
  	    # (onscreen part of the matrix).
  	    # It is the case when column index is strictly superior to
  	    # row index.
  	    if (coords[2] > coords[1]) {
  	      inter.tooltip <- TRUE
  	      clear.tooltip(inter.tooltip, prime.plot)
  	      draw.tooltip(coords, GGI, leg)
  	    } else {
  	      keep.on <- clear.tooltip(inter.tooltip, prime.plot)
  	      inter.tooltip <- !keep.on
  	    }
  	  } else {
  	    keep.on <- clear.tooltip(inter.tooltip, prime.plot)
  	    inter.tooltip <- !keep.on
  	  }
  	}
  }
}

# Function that computes the graphic window x and y extreme values.
# No real plotting happens in this function (empty window is opened).
plot.setup <- function(GGI, colbar.width, draw.names, threshold) {
  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.numeric(colbar.width)) {
    stop("colbar.width argument should be numeric.")
  } else if (!is.logical(draw.names)) {
    stop("draw.names argument should be TRUE or FALSE.")
  } else if (!is.null(threshold) && is.character(threshold) && threshold != "R") {
    stop("threshold argument can not be any other string than 'R'.")
  }

  # Widths and heights of elements are calculated
  if (is.null(threshold)) {
    colorLegend.height <- nrow(GGI)
    colorLegend.width  <- max(0.5, (ncol(GGI)/2)*colbar.width)
    matCol.padding <- colorLegend.width * 0.5
    colorLegend.space <- 1
  } else {
    colorLegend.height <- 0
    colorLegend.width <- 0.5
    matCol.padding <- colorLegend.width * 0.5
    colorLegend.space <- 0
  }

  plot.width  <- ncol(GGI) + colorLegend.space + colorLegend.width + matCol.padding
  plot.height <- nrow(GGI)

  plot(0, xlim=c(2, plot.width), ylim=c(1, plot.height), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")

  # If names are to be plotted and exist, text padding is calculated
  if (draw.names & ncol(GGI) <= 25 & !is.null(colnames(GGI))){
    names.length <- strwidth(colnames(GGI))

    text.vpadding <- ceiling(max(sin(pi/4) * names.length[-ncol(GGI)])) + 0.25
    text.lpadding <- floor(min(seq(2, nrow(GGI)) - 0.25 - names.length[-1]))
    text.rpadding <- ceiling(max(cos(pi/4) * names.length[-ncol(GGI)]))
  } else {
    text.vpadding <- 2
    text.lpadding <- 0
    text.rpadding <- 0
  }

  if (is.null(threshold)) {
    colbar.text.padding <- ceiling(colorLegend.width*0.1 + strwidth("0.75"))
  } else {
    colbar.text.padding <- 0
  }

  xlim <- c(text.lpadding, plot.width + colbar.text.padding + text.rpadding)
  ylim <- c(1, plot.height + text.vpadding)
  plot(0, xlim=xlim, ylim=ylim, type="n",
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
}

# Function that draw the upper triangle of GGI matrix.
# Cells are colored according to corresponding p-values and
# the value of threshold.
# Invisibly return the coordinates of the bottom left point of each
# square. (to save some computing time later)
draw.matrix <- function(GGI, col.FUN, threshold, NA.col, R.thresh, use.log = FALSE) {
  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.function(col.FUN)) {
    stop("col.FUN argument should be a function resulting from colorRampPalette.")
  } else if (!is.logical(use.log)) {
    stop("use.log argument should be a logical.")
  } else if (!is.numeric(R.thresh) || length(R.thresh) != 4) {
    stop("R.thresh argument should be a numeric vector of length 4.")
  } else if (!is.character(NA.col)) {
    stop("NA.col argument should be a character.")
  }

  rect.data <- GGI[upper.tri(GGI)]

  # Assigning colors depending on display options
  if (is.null(threshold)){
    # If gradient is displayed then probs are turned into a percentage first
    if (use.log) {
      quantiles <- c(0, max(GGI, na.rm=TRUE))
    } else {
      quantiles <- c(0, 1)
    }

    rect.perc <- (rect.data - quantiles[1]) / (diff(quantiles))

  } else if (is.numeric(threshold)) {
    # If a threshold is used
    if (use.log) {
      rect.perc <- ifelse(rect.data >= threshold, 0, 1)
    } else {
      rect.perc <- ifelse(rect.data <= threshold, 0, 1)
    }
  } else if (is.character(threshold)) {
    if (use.log){
      rect.perc <- findInterval(rect.data, rev(R.thresh))

    } else {
      rect.perc <- findInterval(rect.data, R.thresh)
    }

      rect.perc <- rect.perc/4
  }

  rect.perc <- floor(rect.perc*200)
  rect.perc[rect.perc == 0] <- 1
  rect.col <- col.FUN(200)[rect.perc]

  # NA values are also disabled
  rect.col[which(is.na(rect.data))] <- NA.col

  rect.pos  <- which(upper.tri(GGI), arr.ind = TRUE)
  temp.X <- rect.pos[, 2]
  rect.pos[, 2] <- max(rect.pos[, 1]) - rect.pos[, 1] + 1
  rect.pos[, 1] <- temp.X

  rect(xleft = rect.pos[, 1],
       ybottom = rect.pos[, 2],
       xright  = rect.pos[, 1] + 1,
       ytop = rect.pos[, 2] + 1,
       col  = rect.col)

  invisible(rect.pos)
}

# Function that draws the gradient indicator.
# The gradient bar is sliced accross the height into a
# large number of smaller rectangles.
draw.colbar <- function(GGI, col.FUN, colbar.width, threshold, R.thresh, NA.col, use.log = FALSE) {
  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.function(col.FUN)) {
    stop("col.FUN argument should be a function resulting from colorRampPalette.")
  } else if (!is.logical(use.log)) {
    stop("use.log argument should be a logical.")
  } else if (!is.numeric(R.thresh) || length(R.thresh) != 4) {
    stop("R.thresh argument should be a numeric vector of length 4.")
  } else if (!is.character(NA.col)) {
    stop("NA.col argument should be a character.")
  }

  if (is.null(threshold)){
    colorLegend.height <- nrow(GGI)
    colorLegend.width  <- max(0.5, (ncol(GGI)/2)*colbar.width)
    matCol.padding <- colorLegend.width * 0.5

    NA.height <- 0.05 * colorLegend.height
    NA.padding <- 0.5 * matCol.padding
    colbar.start <- NA.height + NA.padding

    rect(xleft = rep(ncol(GGI) + 1 + matCol.padding, 200),
         ybottom = seq(1 + colbar.start, ncol(GGI), length=201)[-201],
         xright = rep(ncol(GGI) + 1 + matCol.padding + colorLegend.width, 200),
         ytop = seq(1 + colbar.start, ncol(GGI), length=201)[-1],
         col = col.FUN(200),
         border = NA)

    rect(xleft = ncol(GGI) + 1 + matCol.padding,
         ybottom = 1 + colbar.start,
         xright = ncol(GGI) + 1 + matCol.padding + colorLegend.width,
         ytop = ncol(GGI),
         col = NA,
         border = "black")

    if (use.log) {
      quantiles <- as.numeric(format(quantile(GGI, na.rm=TRUE, names=FALSE), digits=2))
      quantiles.pos <- seq(0, 1, 0.25)
    } else {
      quantiles <- c(0, 0.05, 0.25, 0.5, 0.75, 1)
      quantiles.pos <- quantiles
    }

    segments(x0 = rep(ncol(GGI) + 1 + matCol.padding + colorLegend.width, 6),
             y0 = (ncol(GGI) - 1 -colbar.start) * quantiles.pos + 1 + colbar.start,
             x1 = rep(ncol(GGI) + 1 + matCol.padding + colorLegend.width*1.1 , 6),
             y1 = (ncol(GGI) - 1 -colbar.start) * quantiles.pos + 1 + colbar.start)

    text(x = rep(ncol(GGI) + 1 + matCol.padding + colorLegend.width*1.1 , 6),
         y = (ncol(GGI) - 1 - colbar.start) * quantiles.pos + 1 + colbar.start,
         labels = quantiles,
         pos = 4)

    # NA legend
    rect(xleft = ncol(GGI) + 1 + matCol.padding,
         ybottom = 1,
         xright = ncol(GGI) + 1 + matCol.padding + colorLegend.width,
         ytop =  1 +NA.height,
         col = NA.col,
         border = "black")

    text(x = ncol(GGI) + 1 + matCol.padding + colorLegend.width*1.1,
         y = 1 + 0.5*NA.height,
         labels = "NA",
         pos = 4)

    return(NULL)

  } else if (is.numeric(threshold)) {
    sign <- c("<", ">")

    leg <- legend("bottomleft", c(paste(sign, round(threshold, 3)), 'NA'),
           fill = c(col.FUN(200)[c(1, 200)], NA.col)
           )
  } else if (is.character(threshold)) {
    if (!use.log) {
      legends <- c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "> 0.1")
    } else {
      legends <- c("< 1", "> 1", "> 1.3", "> 2", "> 3")
    }

    leg <- legend("bottomleft", c(legends, 'NA'),
           fill = c(col.FUN(200)[c(1, 50, 100, 150, 200)], NA.col)
    )
  }

  return(leg)
}

# Function that draws genes' names on the plot.
# As the diagonale of the matrix is not drawn, the first
# names is skipped vertically and the last horizontally.
draw.genes.names <- function(genes.names, rect.pos) {
  if(!is.list(genes.names) && length(genes.names) != 2) {
    stop("genes.names argument should be a list of length two.")
  } else if (!is.character(genes.names[[1]]) | !is.character(genes.names[[2]])) {
    stop("genes.names argument should be a list of character vectors.")
  } else if (length(genes.names[[1]]) != length(genes.names[[2]])) {
    stop("genes.names[[1]] & genes.names[[2]] should be of same length.")
  } else if (length(genes.names[[1]]) > 25) {
    stop("Can't handle more than 25 names.")
  } else if (!is.matrix(rect.pos) && !is.numeric(rect.pos[1, 1])) {
    stop("rect.pos argument should be a numeric matrix")
  }

  cex=sort((1 - 1/3:25), decreasing=TRUE)[length(genes.names[[1]]) - 2]

  # Horizontaly
  text(x = sort(unique(rect.pos[, 1])) - 0.25,
       y = sort(unique(rect.pos[, 2]), decreasing = TRUE) + 0.5,
       labels = genes.names[[1]][-length(genes.names[[2]])],
       pos = 2)

  # Verticaly
  text(x = sort(unique(rect.pos[, 1])) + 0.5*min(c(1, (1/length(genes.names[[1]]) ))),
       y = max(rect.pos[, 2]) + 1 + 0.25,
       labels = genes.names[[2]][-1],
       pos = 4,
       srt = 45)
}

# Function that draws the interaction p-values on the matrix.
# Multiple cex are tested so that it is ensured that p-values
# fit inside the squares and are still big enough.
draw.interp <- function(GGI, rect.pos){
  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.matrix(rect.pos) && !is.numeric(rect.pos[1, 1])) {
    stop("rect.pos argument should be a numeric matrix")
  }

  rect.data <- GGI[upper.tri(GGI)]
  rect.data <- format(rect.data, digits=2, scientific=TRUE)

  for (i in seq(1, 0, length=30)) {
    cex <- i
    pval.width <- strwidth(rect.data, cex=cex)
    if (max(pval.width) < 0.9) { break }
  }

  x <- rect.pos[, 1] + 0.5
  y <- rect.pos[, 2] + 0.5
  text(x, y, labels=rect.data, cex=cex)
}

# Function that draws the tooltip windows
# A white black-bordered box is first created
# and text is plotted on top of it.
draw.tooltip <- function(coords, GGI, legend.box) {
  if (is.null(legend.box)) {
    bottomleft <- par('usr')[c(1, 3)]
  } else  {
    bottomleft <- c(legend.box$left + legend.box$w, legend.box$top - legend.box$h)
  }

  tooltip.str <- paste0('Interaction: ',
                       rownames(GGI)[coords[1]], ':',
                       colnames(GGI)[coords[2]],
                       '\np-val: ',
                       format(GGI[coords[1], coords[2]],
                              digits=4))

  rect(xleft = bottomleft[1],
       ybottom = bottomleft[2],
       xright = bottomleft[1] + strwidth(tooltip.str)*1.1,
       ytop = bottomleft[2] + strheight(tooltip.str)*1.5,
       col = "white")
  text(x = bottomleft[1] + strwidth(tooltip.str)*0.05,
       y = mean(c(bottomleft[2], bottomleft[2] + strheight(tooltip.str)*1.5)),
       labels = tooltip.str,
       pos = 4, offset=0)

}

# Function that handles tooltip clearing and tooltip
# procedure ending.
clear.tooltip <- function(inter.tip, prime.plot) {
  # If not in upper triangle then tooltip if cleared
  if (inter.tip) {
    replayPlot(prime.plot)
    return(TRUE)
  } else {
    # If tooltip is already cleared then interaction with
    # user is ceased.
    return(FALSE)
  }
}


draw.network <- function(GGI,genes=1:ncol(GGI),threshold=0.05,plot.nointer=TRUE,method.adjust=c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr")){
  if(length(genes)<2 || length(genes)>ncol(GGI)){
    stop("Number of genes selected not valid.")
  } else if(!class(GGI)%in%c("data.frame","matrix")){
    stop("GGI must be a data.frame.")
  } else if(ncol(GGI)!=nrow(GGI)){
    stop("GGI must be a sqared matrix, containing the pValues for each interaction between genes.")
  } else if(!class(threshold)%in%c("numeric","integer","NULL")){
    stop("Threshold must be a numeric.")
  } else if(class(threshold)=="NULL"){
  	threshold<-0.05
  } else if(threshold>1 || threshold<0){
    stop("Threshold must be comprised in [0,1].")
  } else if(class(plot.nointer)!="logical"){
    stop("plot.inter must be a boolean.")
  }

  if(class(genes)=="character"&&any(!genes%in%colnames(GGI))){
    stop("Genes and GGI don't match. Please select genes that are named in GGI.")
  }

  GGI <- GGI[genes,genes]
  
  method.adjust <- match.arg(method.adjust)
  
  GGI[lower.tri(GGI)] <- p.adjust(GGI[lower.tri(GGI)],method=method.adjust)
  GGI[upper.tri(GGI)] <- p.adjust(GGI[upper.tri(GGI)],method=method.adjust)

  dim <- ncol(GGI)
  pVal.raw <- GGI[lower.tri(GGI)]
  if(any(is.na(pVal.raw))){
    warning("NAs found in GGI, considered as not significative.")
    pVal.raw[is.na(pVal.raw)]<-1
  }

  from.raw <- c()
  to.raw <- c()

  for (i in 1:(dim-1)){
    from.raw <- c(from.raw, rep(colnames(GGI)[i], dim-i))
    to.raw <- c(to.raw, rownames(GGI)[(i+1):dim])
  }



  from <- from.raw[pVal.raw<threshold]
  to <- to.raw[pVal.raw<threshold]
  pVal <- pVal.raw[pVal.raw<threshold]

  if(plot.nointer){
    actors <- data.frame(name=levels(as.factor(unique(c(from.raw,to.raw)))))
  } else {
    actors <- data.frame(name=levels(as.factor(unique(c(from,to)))))
  }

  if(length(from)==0){warning("No interactions has been found between the genes selected.")}

  relations <- data.frame(from=from,
                          to=to,
                          pVal=pVal)

  g <- igraph::graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  plot(g, vertex.size=10)

}


