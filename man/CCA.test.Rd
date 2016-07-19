\name{CCA.test}
\alias{CCA.test}
\title{CCA (Canonical-Correlation Analysis) based GGI analysis.}
\usage{
CCA.test(Y, G1, G2, n.boot = 500)
}
\arguments{
\item{Y}{numeric or factor vector with exactly two different values. \code{Y} is
the response variable and should be of length equal to the number of rows of
\code{G1} and \code{G2} arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{n.boot}{positive integer. \code{n.boot} is the number of bootstrap
replicates for estimating variances. By default, this is fixed to 500.}
}
\value{
A list with class \code{"GGItest"} containing the following components:

	\item{statistic}{The value of the statistic CCU}
	\item{p.value}{The p-value for the test}
	\item{method}{The method used to perform the Gene-Gene interaction test}
	\item{parameter}{A list of the parameters used to perform the Gene-Gene Interaction test}
}
\description{
\code{CCA.test} performs a Gene-Gene Interaction (GGI) analysis based on the
difference of canonical correlation between cases and controls.
}
\details{
	The test statistic is based on the difference between Fisher's transformation of the maximum of the canonical correlations in cases and controls. To calculate the test statistic for the interaction pvalue, \code{CCA.test} estimates the variance of the Fisher's transformation of the maximum of the canonical correlations in cases and controls using a bootstrap method.
}
\examples{
data(gene.pair)
CCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2)
}
\references{
Qianqian Peng, Jinghua Zhao, and Fuzhong Xue. A gene-based method for detecting gene-gene co-association in a case-control study. European Journal of Human Genetics, 18(5) :582-587, May 2010.\cr
}
\seealso{
\code{\link{GGI}}, \code{\link{KCCA.test}}
}

