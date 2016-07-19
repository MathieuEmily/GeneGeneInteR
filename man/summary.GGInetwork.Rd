\name{summary.GGInetwork}
\alias{summary.GGInetwork}
\title{Printing summaries class \code{"GGInetwork"} objects}
\usage{
\method{summary}{GGInetwork}(object, \dots)
}
\arguments{
\item{object}{\code{GGInetwork}. Given as output of \code{\link{GGI}}.}
\item{\dots}{further arguments passed to or from other methods}
}
\description{
	Printing summaries for objects of class \code{"GGInetwork"}
}
\details{
\code{summary.GGItest} provides a customized summary of a \code{GGInetwork} object.
}
\value{
The form of the value returned by \code{summary} depends on the class of its argument. 
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

summary(GGI.res)

}

