\name{imputeSnpMatrix}
\alias{imputeSnpMatrix}
\title{Missing values handling for SnpMatrix object}
\usage{
imputeSnpMatrix(snpX, genes.info, on.rem = c("SNP", "ind", "none"),
  quiet = FALSE)
}
\arguments{
\item{snpX}{\code{snpMatrix} object of which SNPs are to be removed}

\item{genes.info}{A \code{\link[base]{data.frame}} with four columns named \code{Genenames},
\code{SNPnames}, \code{Position} and \code{Chromosome}. Each row describes a
SNP and missing values are not allowed.}

\item{on.rem}{\emph{(optional)} a character string matching one of the
following items: SNP, ind, none. Describes the action taken in case of
remaining missing values. \emph{See details}}

\item{quiet}{\emph{(optional)} a boolean describing wether or not progress bar
should be displayed.}
}
\value{
a list object with two named elements:
\item{\code{snpX}}{a \code{SnpMatrix} object corresponding to input matrix with imputed values.}
\item{\code{genes.info}}{a data frame object corresponding to an updated
 version of input \code{genes.info} in case SNP had to be removed.}

 A warning is printed when SNPs or individuals are removed.
}
\description{
\code{imputeSnpMatrix} is a generic wrapper of
 \code{snp.imputation} and \code{impute.snps} functions from \code{snpStats}
 package. This function mimics a Leave-One-Out process where missing SNPs are
 imputed for an individual based on a model trained on all other individuals.
}
\details{
For the ith row in the \code{snpX} argument (i.e. the ith individual in the dataset), the following steps are performed:
\itemize{
\item missing SNPs are detected for individual i
\item rules are imputed for every missing SNPs using the whole dataset where individual i is removed
\item SNPs are imputed for individual i
}
Although, such a process allows the imputation of a large number of missing SNPs, some missing values may remains. In that case, the user can specify the action to be done thanks to the \code{om.rem} arguments: 
\itemize{
\item \code{om.rem="none"}: leave the dataset as it is, 
\item \code{om.rem="SNP"}: remove all SNPs with remaining missing values,
\item \code{om.rem="ind"}: remove all individuals with remaining missing values.
}
Removing all SNPs is often more parsimonious than removing individuals and allows to get a dataset without any missing values with minimum information-loss.
}
\examples{
\dontrun{
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")
data.imported <- importFile(file=ped, snps=info, pos=posi, pos.sep="\t")

## Example 1 without SNPs filtering
## In this example, 9 SNPs were removed due to remaining missing values.
imputed.snps.1 <- imputeSnpMatrix(data.imported$snpX, genes.info = data.imported$genes.info)


## Example 2 with SNPs filtering priori to the imputation
## Filtering of the data: SNPs with MAF < 0.05 or p.value for HWE < 1e-3 or SNPs with
## call.rate < 0.9 are removed. 
data.scour <- snpMatrixScour(snpX=dta$snpX,genes.info=dta$genes.info,min.maf=0.05,
                              min.eq=1e-3,call.rate=0)

}

## Equivalent loading of data for example 2
## Imputation of the missing genotypes
load(system.file("extdata/dataScour.Rdata", package="GeneGeneInteR"))

data.imputed <- imputeSnpMatrix(data.scour$snpX, genes.info = data.scour$genes.info)
}
\seealso{
\code{\link{GGI}} \code{\link{snpMatrixScour}}
}

