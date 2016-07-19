\name{summary.GGItest}
\alias{summary.GGItest}
\title{Printing summaries class \code{"GGItest"} objects}
\usage{
\method{summary}{GGItest}(object, \dots)
}
\arguments{
\item{object}{\code{GGItest}. Given as output of \code{\link{CCA.test}} for example.}
\item{\dots}{further arguments passed to or from other methods}
}
\description{
	Printing summaries for objects of class \code{"GGItest"}
}
\details{
\code{summary.GGItest} provides a customized summary of a \code{GGItest} object.
}
\value{
The form of the value returned by \code{summary} depends on the class of its argument. 
}
\examples{

data(gene.pair)
summary(PCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,threshold=0.7,method="Std"))


}

