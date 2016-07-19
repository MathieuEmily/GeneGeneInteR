#' Diallelic SNPs informations.
#'
#' A dataset containing datas about 429 individuals and their genotypes.
#'
#' @format A list with 3 objects: \describe{
#'   \item{snpX}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix}
#'    object with 312 SNPs and 429 individuals.} \item{genes.info}{a data frame where
#'   each SNP is described by its rs ID, its position and the gene it belongs
#'   to.} \item{Y}{Factor vector of length 429 with 2 levels (Control or
#'   Rheumatoid - case-), response variable.}}
#' @source All three objects were taken from NCBI web site and are part of the
#'   \href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39428}{GSE39428
#'   series}
"data.SNP"
