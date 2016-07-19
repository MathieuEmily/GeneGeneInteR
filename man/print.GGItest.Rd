\name{print.GGItest}
\alias{print.GGItest}
\title{Printing objects of class \code{"GGItest"}}
\usage{
\method{print}{GGItest}(x, \dots)
}
\arguments{
\item{x}{\code{GGItest}. Given as output of \code{\link{CCA.test}} for example.}
\item{\dots}{further arguments passed to or from other methods} 
}
\value{
the argument x, invisibly, as for all \code{\link{print}} methods.}
\description{
	Printing objects of class \code{"GGItest"}
}
\details{
\code{print.GGItest} provides a customized output of a \code{GGItest} object.
}
\examples{

data(gene.pair)
print(PCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,threshold=0.7,method="Std"))

}

