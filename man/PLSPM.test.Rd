\name{PLSPM.test}
\alias{PLSPM.test}
\title{PLSPM (Partial Least Squares Path Modelling Analysis) based Gene-Gene Interaction analysis.}
\usage{
PLSPM.test(Y, G1, G2,n.perm=500)
}
\arguments{
\item{Y}{numeric or factor vector with exactly two different values. \code{Y} is
the response variable and should be of length equal to the number of rows of
\code{G1} and \code{G2} arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{n.perm}{positive integer. \code{n.perm} is the number of permutations performed to compute the pvalue. By default, this is fixed to 500.}
}
\value{
A list with class \code{"htest"} containing the following components:

	\item{statistic}{The value of the statistic U.}
	\item{p.value}{The p-value for the test.}
	\item{estimate}{A vector of the path coefficients in cases and controls.}
	\item{parameter}{The number of boostrap samples used to estimate the p-value.}
	\item{null.value}{The value of U under the null hypothesis.}
	\item{alternative}{a character string describing the alternative.}
	\item{method}{a character string indicating the method used.}
	\item{data.name}{a character string giving the names of the data.}	
}
\description{
\code{PLSPM.test} performs a Gene-Gene Interaction (GGI) analysis based on the
modelisation of a netwrok of statistical relations. The aim is to quantify the
connections between the latent and the manifest variables.
}
\details{
	The PLSPM-based method, as described in Zhang et al. (2013), aims at comparing the path coefficients \eqn{\beta_D} and \eqn{\beta_C}, where \eqn{\beta_D} is calculated among cases and \eqn{\beta_C} among controls. The co-association between genes \code{G1} and \code{G2} is defined by: \deqn{U=\frac{\beta_D-\beta_C}{\sqrt{Var(\beta_D-\beta_C)}}.}
	The \code{\link[plspm]{plspm}} R package is used to estimate \eqn{U}. The significance pvalue is obtained by using a permutation method on the difference between the path coefficients.
}
\examples{
data(gene.pair)
PLSPM.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,n.perm=50)
}
\references{
X. Zhang et al. (2013) A PLSPM-based test statistic for detecting gene-gene co-association in genome-wide association study with case-control design. PLoS ONE, 8(4) :e62129.\cr
}
\seealso{
\code{\link{GGI}}
}

