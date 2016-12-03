\name{CLD.test}
\alias{CLD.test}
\title{CLD (Composite Linkage Desequilibrium) based GGI analysis.}
\usage{
CLD.test(Y, G1, G2, n.perm = 1000)
}
\arguments{
\item{Y}{numeric or factor vector with exactly two different values. \code{Y} is
the response variable and should be of length equal to the number of rows of
\code{G1} and \code{G2} arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{n.perm}{positive integer. \code{n.perm} is the number of permutations performed to compute the pvalue. By default, \code{n.perm} is fixed to 1000.}
}
\value{
A list with class \code{"htest"} containing the following components:

	\item{statistic}{The value of the statistic CLD.}
	\item{p.value}{The p-value for the test.}
	\item{estimate}{The estimation of CLD}
	\item{parameter}{The number of permutations used to estimate the p-value.}
	\item{null.value}{The value of CLD under the null hypothesis.}
	\item{alternative}{a character string describing the alternative.}
	\item{method}{a character string indicating the method used.}
	\item{data.name}{a character string giving the names of the data.}	
}
\description{
\code{CLD.test} performs a Gene-Gene Interaction (GGI) analysis based on the
difference between the Composite Linkage Disequilibrium measured in cases and controls respectively.
}
\details{
The test statistic is based on Nagao normalized
Quadratic Distance (NQD), as described in Rajapakse et al. (2012). The pvalue is calculated using \code{n.perm} permutations of \code{Y}.
}
\examples{
data(gene.pair)
CLD.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2)
}
\references{
Indika Rajapakse, Michael D. Perlman, Paul J. Martin, John A. Hansen, and Charles Kooperberg. Multivariate detection of gene-gene interactions. Genetic Epidemiology, 36(6) :622-630, 2012.\cr
}
\seealso{
\code{\link{GGI}}
}

