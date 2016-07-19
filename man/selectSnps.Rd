\name{selectSnps}
\alias{selectSnps}
\title{Selects SNPs in snpMatrix and associated info data.frame, outputs of importFile.}
\usage{
selectSnps(snpX, genes.info, select)
}
\arguments{
\item{snpX}{\href{http://bioconductor.org/packages/snpStats/}{snpMatrix} object. Given as output of \code{\link{importFile}}.}

\item{genes.info}{Data.frame containing informations about snps. For more
details, refer to \code{\link{importFile}} help file.}

\item{select}{Numeric or character vector for selecting snps in \code{snpX}
and \code{genes.info}. See details for more information.}
}
\value{
A list of two objects : 

\item{snpX}{a \href{http://bioconductor.org/packages/snpStats/}{snpMatrix}}
\item{genes.info}{ a data frame with 4 columns, and one row per SNP selected with \code{select}. The columns are Chromosome, Genenames, SNPnames and Position.}

An error message is displayed if the genes of snps selected are not found in the either \code{snpX} or \code{genes.info}.
}
\description{
\code{selectSnps} allows the user to select snps from an object output of
\code{\link{importFile}}. It generates the same object, with the columns of
the \href{http://bioconductor.org/packages/snpStats/}{snpMatrix} and the rows
of the data.frame corresponding to the selected snps.
}
\details{
The column names of the \code{genes.info} data.frame should correspond to the output \code{genes.info} object returned by \code{\link{importFile}} function.

The \code{select} argument should one of the following: 
\itemize{
\item a numeric vector with only the column number in the \href{http://bioconductor.org/packages/snpStats/}{snpMatrix} (or row number
for \code{genes.info}) of each snp selected.
\item a character vector with the names of each snp selected or each gene selected.
\item a character vector which elements are position bounds of genes. Each element of the vector is either of the form "begin:end", or "chr:begin:end" if you have to precise the chromosome of the gene.}
}
\examples{
## Importation of the dataset
\dontrun{
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")
data.imported <- importFile(file=ped, snps=info, pos=posi, pos.sep="\t")
}
### Equivalent loading of the imported data
load(system.file("extdata/dataImported.Rdata", package="GeneGeneInteR"))

## Selection of the genes DNAH9 and TXNDC5
selec <- selectSnps(data.imported$snpX, data.imported$genes.info, c("DNAH9","TXNDC5"))

## Selection of the snps from position 101342000 to 101490000 on chromosome 15
selec <- selectSnps(data.imported$snpX, data.imported$genes.info, c("15:101342000:101490000"))
}

