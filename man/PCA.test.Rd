\name{PCA.test}
\alias{PCA.test}
\title{PCA-based Gene-Gene Interaction analysis.}
\usage{
PCA.test(Y, G1, G2, threshold = 0.8,method="GenFreq")
}
\arguments{
\item{Y}{numeric or factor vector with exactly two different values. \code{Y} is
the response variable and should be of length equal to the number of rows of
\code{G1} and \code{G2} arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of rows equal to the length of \code{Y}.}

\item{threshold}{(optional) numeric comprised in \eqn{[0, 1]} interval. \code{threshold} is the minimum percentage of inertia kept for each PCA.}

\item{method}{(optional) character string for PCA method. Must be one of the following: "GenFreq", "Std" (See details).}
}
\value{
A list with class \code{"htest"} containing the following components:

	\item{statistic}{The value of the statistic: the deviance of the anova test.}
	\item{p.value}{The p-value for the test.}
	\item{estimate}{A vector of the residual deviances.}
	\item{parameter}{The degrees of freedom of the chi-squared distribution of the test statistic.}
	\item{null.value}{the value of the deviance under the null.}
	\item{alternative}{a character string describing the alternative.}
	\item{method}{a character string indicating the method used.}
	\item{data.name}{a character string giving the names of the data.}	
}
\description{
\code{PCA.test} performs a Gene-Gene Interaction (GGI) analysis by testing the interaction between the principal components of the two genes. With \code{method="Std"} PCA is standardized using standard deviation for each variable. With \code{method="GenFreq"}, dataset is standardized using the standard deviation under Hardy-Weinberg equilibrium, as proposed in the \href{http://bioconductor.org/packages/snpStats/}{snpStats} Bioconductor package. 
}
\details{
In a first step, \code{PCA.test} performs a Principal Components Analysis on both \code{G1} and \code{G2} genes that are interpreted as matrices of allele counts. With the \code{method="Std"}, the dataset is standardized using variables standard deviation,  while dataset is standardized using standard deviation under Hardy-Weinberg equilibrium for \code{method="GenFreq"}.  Principal components are then retrieved to describe each dataset with user-defined inertia percentage (parameter \code{threshold}) and used in a logistic regression model. The consists is testing the significance of the interaction terms using a Likelihood Ratio Test (see Li et al. (2009)).}
\examples{
data(gene.pair)
PCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,threshold=0.7,method="Std")
PCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,threshold=0.7,method="GenFreq")
}
\references{
J. Li et al. (2009) Identification of gene-gene interaction using principal components. BMC Proceedings, 3 (Suppl. 7): S78\cr
}
\seealso{
\code{\link{GGI}}
}

