\docType{data}
\name{data.SNP}
\alias{data.SNP}
\title{Multiple genes dataset}
\format{A list with 3 objects: \describe{
  \item{snpX}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix}
   object with 312 SNPs and 429 individuals.} \item{genes.info}{a data frame where
  each SNP is described by its rs ID, its position and the gene it belongs
  to.} \item{Y}{Factor vector of length 429 with 2 levels ("Health Control or
  "Rheumatoid"). \code{Y} correspond to the case-control status response variable.}}}
\source{
All three objects were taken from NCBI web site and are part of the
  \href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39428}{GSE39428
  series}
Chang X., et al. Investigating a pathogenic role for TXNDC5 in tumors. Int. J. Oncol.,  43(6): 1871-84, 2013.\cr
}
\usage{
data(data.SNP)
}
\description{
A case-control dataset containing the genotypes of 312 SNPs from 17 genes in a total of 429 patients (266 individuals affected by Rheumatoid Arthritis and 163 Health controls)
}
\value{
A dataset containing the genotypes of 312 SNPs from 17 genes in a total of 429 patients.
}
\keyword{datasets}

