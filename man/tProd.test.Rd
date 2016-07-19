\name{tProd.test}
\alias{tProd.test}
\title{Gene-based Gene-Gene Interaction analysis by combining SNP-SNP interaction tests}
\usage{
tProd.test(Y, G1, G2, tau = 0.05, n.sim = 1000)
}
\arguments{
\item{Y}{numeric or factor vector with exactly two different values. \code{Y} is
the response variable and should be of length equal to the number of rows of
\code{G1} and \code{G2} arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{tau}{numeric in [0, 1]. See details section for its use.}

\item{n.sim}{positive integer. \code{n.sim} is the number of multivariate normal distribution simulations used to compute the significant level for tTS and tProd methods.}
}
\value{
A list with class \code{"GGItest"} containing the following components:

	\item{statistic}{The value of the statistic tTS}
	\item{p.value}{The p-value for the test}
	\item{method}{The method used to perform the Gene-Gene interaction test}
	\item{parameter}{A list of the parameters used to perform the Gene-Gene Interaction test}

}
\description{
\code{gates.test}, \code{minP.test}, \code{tTS.test} and \code{tProd.test} aim at performing gene-gene interaction analysis based on SNP-SNP interaction tests. The following methods are used to combine SNP-SNP interaction tests
 into a single Gene-Gene Interaction p-value: 
 \itemize{ 
 	\item Minimum p-value in \code{minP.test} function 
 	\item Gene Association Test with Extended Simes procedure in \code{gates.test}
 	\item Truncated Tail Strength in \code{tTS.test} function
 	\item Truncated p-value Product in \code{tProd.test} function
}
}
\details{
	In a first step, all methods start by applying a logistic regression model to test all pairs of SNPs between the two genes \code{G1} and \code{G2}. If \code{G1} has \eqn{m1} SNPs and \code{G2} \eqn{m2} SNPs, a total of \eqn{m1 * m2} SNP-SNP tests are performed. In a second step, the \eqn{m1 * m2} SNP-SNP tests are combined according to their covariance matrix \eqn{\Sigma}. \eqn{\Sigma} is computed as described in the method developped in Emily (2016). The covariance \eqn{\Sigma} is used in each method as follows:

 \itemize{ 
 	\item minP test - minP test considered the significant of the observed minimum p-value. Significance is computed by integrating the multivariate normal distribution with covariance \eqn{\Sigma} as proposed in Conneelly and Boehnke (2008).
 	\item GATES test - The p-value for GATES is the minimum p-value obtained after a multiple testing correction of the SNP-SNP interaction p-values. Correction for multiple testing is defined as \deqn{me * p[i]/me[i]} where \eqn{me} is the effective number of independant tests, \eqn{p[i]} is the i-th top p-values and \eqn{me[i]} is the effective number of
 independant test among the top i p-values. Many methods exist to estimate \eqn{me} and \eqn{me[i]} terms: 
 \itemize{
 	\item Cheverud-Nyholt method (Cheverud, 2001 and Nyholt, 2004)
 	\item Keff method (Moskovina and Schmidt, 2008)
 	\item Li & Ji method (Li and Ji, 2005)
 	\item Galwey method (Galwey, 2009)
 } Details of each method can be found in the references. 
 	\item tTS test - tTS test does not consider only the strongest signal but all signals that are inferior to a given threshold \eqn{\tau}. For these
 p-values, the weighted sum of \deqn{tTS=\sum (1-p[i]*(m1*m2+1)/i)} is computed and
 represents the test statistic. The p-value is calculated using an empirical distribution of \eqn{tTS} obtained by simulating multivariate normal statistics with a covariance \eqn{\Sigma} as proposed by Jiang et al. (2011). 
 \item TProd test - The procedure is similar to \eqn{tTS} with \eqn{tProd} defined as \deqn{tProd=\prod p[i].}
See Zaykin et al. (2002) for details.}
}
\examples{
data(gene.pair)
tProd.test(gene.pair$Y, gene.pair$G1, gene.pair$G2, tau = 0.5, n.sim = 500)
}
\references{
M. Emily AGGrEGATOr: A Gene-based GEne-Gene interActTiOn test for case-control association studies, Statistical Application in Genetics and Molecular Biology, 15(2): 151-171, 2016.\cr
L. Ma, A.G. Clark and A. Keinan Gene-Based Testing Of Interactions in Association Studies of Quantitative Traits. PLoS Genetics 9(2):e1003321, 2013.\cr
V. Moskvina and K.M. Schmidt On multiple-testing correction in genome-wide association studies. Genetic Epidemiology, 32(6): 567-573, 2008.\cr
J. Li and L. Ji. Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. Heredity 95: 221-227, 2005.\cr
N.W. Galwey. A new measure of the effective number of tests, a practical tool for comparing families of non-independent significance tests. Genetic Epidemiology 33(7): 559-568, 2009.\cr
J.M. Cheverud. A simple correction for multiple comparisons in interval mapping genome scans. Heredity. 87(1):52-58, 2001.\cr
D.R. Nyholt. A Simple Correction for Multiple Testing for Single-Nucleotide Polymorphisms in Linkage Disequilibrium with Each Other. American journal of human genetics. 74(4): 765-769, 2004.\cr
K.N. Conneely and M. Boehnke. So many correlated tests, so little time! rapid adjustment of p values for multiple correlated tests. The American Journal of Human Genetics, 81: 1158-1168, 2008.\cr
B. Jiang, X. Zhang, Y. Zuo and G. Kang. A powerful truncated tail strength method for testing multiple null hypotheses in one dataset. Journal of Theoretical Biology 277: 67-73, 2011.\cr
D.V. Zaykin, L.A. Zhivotovsky, P.H. Westfall and B.S. Weir. Truncated product method for combining P-values. Genetic epidemiology 22: 170-185, 2002.\cr
}
\seealso{
\code{\link{GGI}}
}

