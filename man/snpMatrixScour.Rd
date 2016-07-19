\name{snpMatrixScour}
\alias{snpMatrixScour}
\title{SNP filtering based on Minor Allele Frequency, Hardy-Weinberg Equilibrium and call rate}
\usage{
snpMatrixScour(snpX, genes.length = NULL, genes.info = NULL,
  min.maf = 0.01, min.eq = 0.01, call.rate = 0.9)
}
\arguments{
\item{snpX}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object from which SNPs are to be removed}

\item{genes.length}{\emph{(optional)} numeric vector. It is the length (in columns/SNPs) of each gene. Each gene declared is considered contiguous with the one before and after it. \code{genes.lengths} can be named (names will be kept).}

\item{genes.info}{\emph{(optional)} a \code{\link[base]{data.frame}} with four columns named \code{Genenames}, \code{SNPnames}, \code{Position} and \code{Chromosome}. Each row describes a SNP and missing values are not allowed.}

\item{min.maf}{a single value between 0 and 0.5 that gives the threshold for the MAF (Minor Allele Frequency) of a SNP. SNP with MAF < \code{min.maf} are removed. Default is 0.01.}

\item{min.eq}{a single value between 0 and 1 that gives the maximum acceptable p-value for the \eqn{\chi ^2} verifying HWE deviation. SNP that does not meet that criterion are removed. Default is 0.01.}

\item{call.rate}{a single value between 0 and 1 that gives the minimum acceptable call rate for a SNP. Default is 0.9. Low values for SNPs call rate can make imputation harder (residual missing values).}
}
\value{
A list with two objects:
\describe{
\item{\code{snpX}}{the \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object where non-conform SNPs are removed.}
\item{\code{genes.info}}{the object that contains the updated gene lengths information. Can be a numeric vector (possibly named) or a data frame. If \code{genes.length} and \code{genes.info} are not provided by the user as input of the \code{snpMatrixScour} function, the \code{genes.info} object is \code{NULL}.}
}
}
\description{
\code{snpMatrixScour} aims at filtering out SNPs of a \href{http://bioconductor.org/packages/snpStats/}{snpMatrix} object based on Minor Allele Frequency criterion, deviation to Hardy-Weinberg Equilibrium and SNPs call rate.
}
\details{
This function removes SNPs that does not meet all following criteria:
\itemize{
	\item \eqn{MAF < }\code{min.maf},
	\item \eqn{HWE < }\code{min.eq}, where HWE is the p-value of the \eqn{\chi ^2} test of deviation to Hardy-Weinberg Equilibrium,
	\item Call rate \eqn{> }\code{call.rate}.
}

If \code{genes.length} and \code{genes.info} are provided by the user, an updated version is returned by \code{snpMatrixScour}. The returned object can be directly used as inputs of the \code{\link{GGI}} function.
}
\examples{
\dontrun{
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")
data.imported <- importFile(file=ped, snps=info, pos=posi, pos.sep="\t")
}
### Equivalent loading of the imported data
load(system.file("extdata/dataImported.Rdata", package="GeneGeneInteR"))

## In this example, SNPs are with MAF lower than 0.2 or p-value for HWE testing lower than 0.05 or
# a proportion of missing value higher than 0.2 are removed
data.scour1 <- snpMatrixScour(data.imported$snpX, genes.info = data.imported$genes.info,
                               min.maf = 0.2, min.eq=0.05, call.rate = 0.8)
## Two genes have been completely removed from the resulting dataset.
}
\seealso{
\code{\link{GGI}}
}

