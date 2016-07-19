\docType{data}
\name{gene.pair}
\alias{gene.pair}
\title{A case-control dataset with two genes}
\format{A list with 3 objects:
\describe{
\item{Y}{A \code{\link[base]{factor}} of disease status with two levels: \code{HealthControl} or \code{RheumatoidArthritis}.}
\item{G1}{A \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object with 8 SNPs.}
\item{G2}{A \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object with 4 SNPs.}
}
}
\source{
All three objects were taken from NCBI web site and are part of the  \href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39428}{GSE39428 series}
Chang X., et al. Investigating a pathogenic role for TXNDC5 in tumors. Int. J. Oncol.,  43(6): 1871-84, 2013.\cr
}
\usage{
data(gene.pair)
}
\description{
A case-control dataset containing the genotypes of 8 SNPs within \emph{GC} gene (object \code{G1}) and 4 SNPs within \emph{PADI2} gene (object \code{G2}). The dataset includes 247 individuals affected by Rheumatoid Arthritis (RA) and 202 individuals not affected by RA.
}
\value{
A dataset containing the genotypes of 8 SNPs within \emph{GC} gene (object \code{G1}) and 4 SNPs within \emph{PADI2} gene (object \code{G2}).
}
\keyword{datasets}

