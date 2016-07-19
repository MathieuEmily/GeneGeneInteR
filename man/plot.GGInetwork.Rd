\name{plot.GGInetwork}
\alias{plot.GGInetwork}
\title{Plot a Gene-Gene Interaction Study with an heatmap-like visualization}
\usage{
\method{plot}{GGInetwork}(x,method=c("heatmap","network"),threshold=NULL, 
col=c("#D6604D","#104E8B"),colbar.width=0.15, title=NULL,hclust.order=FALSE, 
use.log=FALSE,NA.col="#D3D3D3",draw.pvals=(ncol(x$p.value) <= 15), 
draw.names=(ncol(x$p.value)<= 25),interact=FALSE, 
method.adjust=c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr"), 
genes=1:ncol(x$p.value), plot.nointer=TRUE, \dots)
}
\arguments{
\item{x}{\code{GGInetwork}. Given as output of \code{\link{GGI}}.}

\item{method}{Output graph ("heatmap" for heatmap-like, "network" for network-like). Default is \code{heatmap}.}

\item{threshold}{A numeric between 0 and 1. All p-value strictly greater than
that \code{threshold} are distinctly colored (See \code{NA.col}). Default is \code{NULL} corresponding to no threshold for applied for an heatmap graph and a 0.05 threshold applied for the network representation.}

\item{method.adjust}{correction method for multiple testing as proposed in the \code{\link[stats]{p.adjust}} function. By default, no correction is applied to the p-values.}

\item{genes}{Numeric vector allowing a selection of the genes that will be included in the relations. Default is set to all genes.}

\item{col}{(\emph{Only for \code{method=heatmap}}.) String vector. Marker colors to be used for the gradient. The first
element of the vector is the value for 0 and the last is for 1. If only one
value is parsed, that color is used for 0 and white is automatically used
for 1. Any value compatible with \code{\link{colorRampPalette}} function.}

\item{colbar.width}{(\emph{Only for \code{method=heatmap}}.) A positive number describing the gradient bar width ratio.
That number is used to keep the gradient bar's width steady as the size of
\code{GGI} increases.}

\item{title}{(\emph{Only for \code{method=heatmap}}.) A string used as the plot title. If left as NULL, a generic name
is generated.}

\item{hclust.order}{(\emph{Only for \code{method=heatmap}}.) A boolean. Should a hierachical clustering procedure be performed on \code{GGI} to order the matrix ?}

\item{use.log}{(\emph{Only for \code{method=heatmap}}.) A boolean. Should the p-values be reported in -log10? FALSE by default.}

\item{NA.col}{(\emph{Only for \code{method=heatmap}}.) A string. The color used when a p-value is strictly greater than \code{threshold}.}

\item{draw.pvals}{(\emph{Only for \code{method=heatmap}}.) A boolean. Should p-values be plotted ? Disabled when the size of \code{GGI} exceeds \eqn{15*15}.}

\item{draw.names}{(\emph{Only for \code{method=heatmap}}.) A boolean. Should genes' names be plotted on matrix margins? Disabled when the size of \code{GGI} exceeds \eqn{25*25}.}

\item{interact}{(\emph{Only for \code{method=heatmap}}.) A boolean. Should the plot be clickable ? (See Details for more information). Disabled when open R session is not interactive.}

\item{plot.nointer}{(\emph{Only for \code{method=network}}.) A boolean. Set TRUE if the genes with no interaction should be plotted.}

\item{\dots}{further arguments passed to or from other methods.}
}
\description{
\code{GGI.plot} is a graphical function that allow a heatmap-like or a network-like visualization of a Gene-Gene Interaction study based on a set of genes.
}
\details{
	If \code{method=heatmap}, this function draw the upper half of a Gene-Gene Interaction results matrix without its diagonal. A gradient is created from 0 to 1 (by default from crimson to white) and the matrix cells are colored according to the corresponding p-value.

 By  default, p-values and names are drawn to make matrix reading easier, but
 in case parameter \code{GGI} is large, p-values (and eventually gene names as
 \code{GGI} grows bigger) are not drawn. In that case, the default behavior
 of the function is to start an interactive process where user can click on a
 cell of interest to open a tooltip displaying which genes are involved in
 selected interaction and the p-value of the interaction test. Tooltips can
 be closed if user clicks anywhere else than on a cell. This process stops
 when the user presses the escape button (or terminates the locator procedure
 in general) or when the user clicks on any place other than a cell when no
 tooltip window is open.

 To improve plot clarity, user may set a threshold above which cells are
 colored with a distinct color. By default, threshold is set to 1 and no
 cell is colored differently (as values must be \emph{strictly} above the
 threshold).
 
 If \code{method=network}, this function plots a graph representing the significative interactions between genes of a Gene-Gene Interaction study.
}
\value{
The form of the value returned by \code{plot} depends on the class of its argument. See Details.
}
\examples{
\dontrun{
## Dataset is included in the package
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")

## Importation of the genotypes
data.imported <- importFile(file=ped, snps=info, pos=posi, pos.sep="\\t")
## Filtering of the data: SNPs with MAF < 0.05 or p.value for HWE < 1e-3 or 
## call rate < 0.9 are removed. 
data.scour <- snpMatrixScour(snpX=data.imported$snpX,genes.info=data.imported$genes.info,
                              min.maf=0.05,min.eq=1e-3,call.rate=0.9)
## Imputation of the missing genotypes
data.imputed <- imputeSnpMatrix(data.scour$snpX, genes.info = data.scour$genes.info)

## Importation of the phenotype
resp <- system.file("extdata/response.txt", package="GeneGeneInteR")
Y  <- read.csv(resp, header=FALSE)

## plot of the interaction between the 17 genes with the CLD method -- can take a few minutes
GGI.res <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info,method="CLD")
plot(GGI.res,threshold=0.05)

## Selection of 12 genes among 17
data.select <- selectSnps(data.imputed$snpX, data.imputed$genes.info, c("bub3","CDSN","Gc","GLRX",
                          "PADI1","PADI2","PADI4","PADI6","PRKD3","PSORS1C1","SERPINA1","SORBS1"))
GGI.res <- GGI(Y=Y, snpX=data.select$snpX, genes.info=data.select$genes.info,method="PCA")
}
## Equivalent importation of the GGI.res object
load(system.file("extdata/GGIRes.Rdata", package="GeneGeneInteR"))


## Plot of the results with default values
plot(GGI.res)

## Plot of the results with a threshold and an ordering of the genes. 
## Default method is an heatmap-like representation
plot(GGI.res,threshold=0.1,hclust.order=TRUE)

## Example of network with default threshold 0.05
plot(GGI.res,method="network")

## Example of network with threshold 0.01 where genes with no interaction are not plotted 
# (plot.nointer=FALSE)
plot(GGI.res,method="network",threshold=0.1,plot.nointer=FALSE)

}
\seealso{
\code{\link{GGI}}
}

