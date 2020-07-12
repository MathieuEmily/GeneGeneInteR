\name{importFile}
\alias{importFile}
\title{Imports SNPs information from pedfile, PLINK, VCF (4.0) file, or genotypes
imputed by IMPUTE2.}
\usage{
importFile(file, pos, pos.sep = "\t", ...)
}
\arguments{
\item{file}{String containing the path of the file to import.}

\item{pos}{(optional) Path of a csv file, character vector or numeric vector
containing informations for each SNP about chromosome, gene names, snp names
and positions.}

\item{pos.sep}{(optional) String to be passed to \code{read.csv} function as
\code{sep} argument. Default is tab (\\t).}

\item{...}{Additionnal arguments to be passed to the reading file function.}
}
\value{
A list of two objects : 
\item{\code{snpX}}{a \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix}}
\item{\code{genes.info}}{a data frame with 4 columns, and one row per SNP. The columns are Chromosome, Genenames, SNPnames and position.}
}
\description{
\code{importFile} generates a
\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object on the
basis of diallelic object contained in a file, and creates also a data frame
containing information about positions of the SNPs on the genome.
}
\details{
As input information, \code{importFile} takes the full paths of the files to be imported. Files are then read with
\code{\link[snpStats]{read.pedfile}}, \code{\link[snpStats]{read.plink}},
\code{\link[GGtools]{vcf2sm}}, or \code{\link[snpStats]{read.impute}},
depending on the file extension (pedfile, plink, vcf or impute2). For a
pedfile, \code{importFile} also reads the ".info" file associated (this must
be in the same directory as the ".ped". Similarly, for a PLINK, there must be
3 files with extensions ".bed" (passed as argument \code{file}), ".bim" and
".fam". A VCF file must be with the associated ".tbi" file.

If the file is a vcf file, two additional arguments have to done : \itemize{\item
\code{gr} instance of \code{\link[GenomicRanges]{GRanges}.} \item
\code{nmetacol} numeric defining the number of columns used in each record as
locus-level metadata.}

\code{Pos} argument is optional. If it's not given, the data frame with position
information is filled with NAs, except for SNP names which are imported from
the column names of the
\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix}, and eventually
positions and chromosomes if the file imported is a pedfile or a PLINK. Else,
the \code{pos} argument can be either the path to a csv file, a character vector with
elements of the form chr:position, or a numeric vector with only the
positions. Additionnaly, SNP names can be precised as names of the vector. If
you choose the csv file path, be sure that the columns are named as
follows : Chromosome, Genenames, SNPnames, Position.}
\examples{
## Pedfile from this package.
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")

## Information about position of the snps
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")

## Importation
data.imported <- importFile(file=ped, snps=info, pos=posi, pos.sep="\t")

#########
## VCF file from GGtools package.
\dontrun{
vref <- system.file("vcf/CEU.exon.2010_09.genotypes.vcf.gz", package="GGtools")
irange <- IRanges::IRanges(10e6,20e6)
gg = GenomicRanges::GRanges(seqnames="1", ranges=irange)
dta <- importFile(file=vref, gr=gg, nmetacol=9L)
}
}

