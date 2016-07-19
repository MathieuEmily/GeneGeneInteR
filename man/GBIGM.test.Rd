\name{GBIGM.test}
\alias{GBIGM.test}
\title{GBIGM (Gene-based information gain method) for GGI analysis.}
\usage{
GBIGM.test(Y, G1, G2, n.perm = 1000)
}
\arguments{
\item{Y}{numeric or factor vector with exactly two different values. \code{Y} is
the response variable and should be of length equal to the number of rows of
\code{G1} and \code{G2} arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{n.perm}{positive integer. \code{n.perm} is the number of permutations performed to compute the pvalue. By default, this is fixed to 1000.}
}
\value{
A list with class \code{"GGItest"} containing the following components:

	\item{statistic}{The value of the statistic \eqn{DeltaR1,2}}
	\item{p.value}{The p-value for the test}
	\item{method}{The method used to perform the Gene-Gene interaction test}
	\item{parameter}{A list of the parameters used to perform the Gene-Gene Interaction test}

}
\description{
\code{GBIGM.test} performs a Gene-Gene Interaction (GGI) analysis by contrasting the information entropy between cases and controls.
}
\details{
First, the conditional entropy and information gain rate are computed for each gene \code{G1} and \code{G2}. In a second step, information gain rate for the gene pair (\code{G1},\code{G2}) is computed. A p-value is estimated using permutations of \code{Y}. More details can be found in Li et al. (2015).
}
\examples{
data(gene.pair)
GBIGM.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,n.perm=500)
}
\references{
J. Li, et al.. A gene-based information gain method for detecting gene-gene interactions in case-control studies. European Journal of Human Genetics, 23 :1566-1572, 2015.\cr
 }

\seealso{
\code{\link{GGI}}
}

