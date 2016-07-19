\name{KCCA.test}
\alias{KCCA.test}
\title{KCCA (Kernel Canonical-Correlation Analysis) based GGI analysis.}
\usage{
KCCA.test(Y, G1, G2,kernel=c("rbfdot","polydot","tanhdot","vanilladot","laplacedot",
"besseldot","anovadot","splinedot"),n.boot = 500,sigma=0.05,degree=1,scale=1,offset=1,
order=1)
}
\arguments{
\item{Y}{numeric or factor vector with exactly two different values. \code{Y} is
the response variable and should be of length equal to the number of rows of
\code{G1} and \code{G2} arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{kernel}{A character string matching one of the kernel name in : "rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot". For more details regarding kernel function see \code{\link[kernlab]{dots}} object. Default is a Gaussian RBF kernel}

\item{n.boot}{positive integer. \code{n.boot} is the number of bootstrap replicates for estimating variances. By default, this is fixed to 500.}

\item{sigma}{The inverse kernel width used by the Gaussian the Laplacian (\code{kernel="rbfdot"}), the Bessel (\code{kernel="besseldot"}) and the ANOVA  (\code{kernel="anovadot"}) kernel. Default is \code{sigma=0.05}.}

\item{degree}{The degree of the polynomial (\code{kernel="polydot"}), bessel (\code{kernel="besseldot"}) or ANOVA (\code{kernel="anovadot"}) kernel function. This has to be an positive integer. Default is \code{degree=1}.}

\item{scale}{The scaling parameter of the polynomial (\code{kernel="polydot"}) and tangent (\code{kernel="tanhdot"}) kernel is a convenient way of normalizing patterns without the need to modify the data itself. Default is \code{scale=1}.}

\item{offset}{The offset used in a polynomial (\code{kernel="polydot"}) or hyperbolic tangent (\code{kernel="tanhdot"}) kernel. Default is \code{offset=1}.}

\item{order}{The order of the Bessel function to be used as a kernel (\code{kernel="besseldot"}). Default is \code{order=1}.} 
}
\value{
A list with class \code{"GGItest"} containing the following components:

	\item{statistic}{The value of the statistic KCCU}
	\item{p.value}{The p-value for the test}
	\item{method}{The method used to perform the Gene-Gene interaction test}
	\item{parameter}{A list of the parameters used to perform the Gene-Gene Interaction test}
}
\description{
\code{KCCA.test} performs a Gene-Gene Interaction (GGI) analysis based on the
difference of canonical correlations between cases and controls. The "kernel trick" is applied to the canonical correlation to allow the detection of non-linear co-association.
}
\details{
	The test statistic is based on the difference between a Fisher's transformation of the maximum of the kernelized canonical correlations in cases and controls. To calculate the test statistic for the interaction pvalue, \code{KCCA.test} estimates the variance of the Fisher's transformation of the maximum of the kernelized canonical correlations in cases and controls using a bootstrap method. The computation of \code{\link[kernlab]{kcca}.} can be very long.
}
\examples{
data(gene.pair)
\dontrun{
KCCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2)
}
}
\references{
Yuan, Z. et al. (2012): Detection for gene-gene co-association via kernel canonical correlation analysis, BMC Genetics, 13, 83. \cr
Larson, N. B. et al. (2013): A kernel regression approach to gene-gene interaction detection for case-control studies, Genetic Epidemiology, 37, 695-703. \cr
}
\seealso{
\code{\link{GGI}}, \code{\link{CCA.test}}
}

