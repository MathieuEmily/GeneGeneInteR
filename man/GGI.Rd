\name{GGI}
\alias{GGI}
\title{Gene-Gene Interaction Analysis of a set of genes}
\usage{
GGI(Y, snpX, genes.length = NULL, genes.info = NULL, method = c("minP","PCA", "CCA",
 "KCCA","CLD","PLSPM","GBIGM","GATES", "tTS", "tProd"), \dots)
}
\arguments{
\item{Y}{numeric or factor vector with exactly two different values. \code{Y} is
the response variable and should be of length equal to the number of rows of
\code{snpX} (number of individuals).}

\item{snpX}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object. Must have a number of rows equal to the length of \code{Y}. See details.}

\item{genes.length}{\emph{(optional)} a numeric vector. \code{gene.length} is the length (in columns/SNP) of each
gene.}

\item{genes.info}{\emph{(optional)} a data frame. \code{genes.info} must have four columns
named \code{Genenames}, \code{SNPnames}, \code{Position} and
\code{Chromosome}. Each row describes a SNP and missing values are not
allowed.}

\item{method}{a string matching one of the following: PCA, CCA, KCCA, CLD, PLSPM, GBIGM, minP, GATES, tTS or tProd. Only one method can be parsed.}

\item{\dots}{Other optional arguments to be passed to the functions associated
with the method chosen. See more in elementary methods help.}
}
\value{
A list with class \code{"GGInetwork"} containing the following components:

\item{statistic}{a symmetric \code{\link{matrix}} of size \eqn{G*G} where \eqn{G} is the number of genes studied. The general term of the matrix is the statistic of the interaction between the two genes.}
\item{p.value}{a symmetric \code{\link{matrix}} of size \eqn{G*G} where \eqn{G} is the number of genes studied. The general term of the matrix is the p-value of the interaction between the two genes.}
\item{df}{\emph{(Only for \code{method="PCA"}).} a symmetric \code{\link{matrix}} of size \eqn{G*G} where \eqn{G} is the number of genes studied. The general term of the matrix is the degrees of freedom of the interaction test.}
\item{method}{The method used to perform the Gene-Gene interaction test.}
\item{parameter}{A list of the parameters used to perform the Gene-Gene Interaction test.}
}
\description{
\code{GGI} allows the search for Gene-Gene Interactions by testing all possible pairs of genes 
in a set of genes.
}
\details{
This function is a wrapper for all Gene-Gene Interaction analysis methods and drive the
 overall analysis: splitting the dataset in gene matrices and starting
 elementary analysis for each pair of genes.

 \emph{SNPs from the same gene are assumed to be ordered along the chromosome. 
 See \code{\link{selectSnps}}.}

 If \code{genes.lenght} is provided, it contains the number of SNPs of each gene. For example, 
 if \code{genes.length} is the vector: c(20, 35, 15), then gene 1 will be interpreted as the set 
 of the first 20 columns/SNPs of \code{snpX}, gene 2 will be interpreted as the following 
 35 columns/SNP, etc. Each gene declared is considered contiguous with the one before and
after it. \code{genes.length} can be named if you want the returned matrix
to have dimensions named after those. If no names are given then generic
names are generated following the pattern Gene.n (n being the gene's index)
.

The following methods are available to perform the interaction test for a single pair of genes:

\itemize{ 
\item Principal Components Analysis method (PCA) \code{\link{PCA.test}} - PCA is performed on both 
genes and resulting principal components are used to fit a logistic regression model with 
interaction between and a second logistic regression model without interaction term. 
The interaction between the two genes is then tested using a likelihood ratio test between the two 
logistic regression models (see Li et al. 2009).
\item Canonical Correlation Analysis (CCA) \code{\link{CCA.test}} - The maximum of canonical 
correlation between the two genes is computed for each group (cases and controls). The difference
 between the two transformed values (Fisher transformation) is used to test for interaction between genes
 (see Peng et al. 2010). 
\item Kernel Canonical Correlation Analysis (KCCA) \code{\link{KCCA.test}} - This method is similar to the CCA method where the canonical correlations are computed using Kernel method (see Yuan et al., 2012
 and Larson et al., 2013). 
\item Composite Linkage Disequilibrium (CLD) \code{\link{CLD.test}} - CLD is  based on the 
difference of the covariance matrices between the two genes computed for cases and controls. 
The covariance is estimated via the Composite Linkage Disequilibrium and a method based on Nagao
 normalized Quadratic Distance is used to compute the test statistic (see Rajapakse et al., 2012). 
\item Partial Least Square Path Modeling (PLSPM) \code{\link{PLSPM.test}} - A network of 
statistical relations between latent and manifest variables is built. The difference between the
 path coefficients is used to compute the test statistic (see Zhang et al., 2013).
\item Gene-Based Information Gain Method (GBIGM) \code{\link{GBIGM.test}} - Entropies and 
Information Gain Ratio are used to compute a measure of the co-association between two genes
 (see Li et al., 2015).
\item Minimum p-value test (minP) \code{\link{minP.test}} - Given two genes, \code{G1} with 
\eqn{m1} SNPs and \code{G2} with \eqn{m2} SNPs, all SNP-SNP interactions are first tested using 
a logistic regression model, thus generated a set of \eqn{m1*m2} p-values. The significance of the
 minimum p-value is evaluated using multivariate normal distribution that accounts for the 
 covariance between the tests statistics at the SNP level (see Emily, 2016).
\item Gene Association Test using Extended Simes procedure (GATES) \code{\link{gates.test}} 
- Given two genes, \code{G1} with \eqn{m1} SNPs and \code{G2} with \eqn{m2} SNPs, all SNP-SNP 
interactions are first tested using a logistic regression model, thus generated a set of
 \eqn{m1*m2} p-values. P-values are then corrected for multiple testing using an extension of the
 Simes procedure that take into account the correlation between the tests statistic via the number
 of effective tests (see Li. et al., 2011). 
\item Truncated Tail Strength test (tTS) \code{\link{tTS.test}} - Given two genes, \code{G1} with
 \eqn{m1} SNPs and \code{G2} with \eqn{m2} SNPs, all SNP-SNP interactions are first tested using 
 a logistic regression model, thus generated a set of \eqn{m1*m2} p-values. All p-values below a 
 user defined threshold are weighted and summed up to provide the tTS test statistic 
 (see Jiang et al., 2011).
\item Truncated p-value Product test (tProd) \code{\link{tProd.test}} - Similar to tTS but with 
a different p-values transformation (see Zaykin, 2002)
}

Missing values are not allowed and trying to parse an incomplete \code{SnpMatrix} object as an 
argument will result in an error. Imputation can be performed prior to the analysis with the
 \code{\link{imputeSnpMatrix}} function.
}
\examples{
\dontrun{
## Dataset is included in the package
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")

## Importation of the genotypes
data.imported <- importFile(file=ped, snps=info, pos=posi, pos.sep="\\t")
## Filtering of the data: SNPs with MAF < 0.05 or p.value for HWE < 1e-3 or SNPs with 
## call.rate < 0.9 are removed. 
data.scour <- snpMatrixScour(snpX=data.imported$snpX,genes.info=data.imported$genes.info,min.maf=0.05,
                              min.eq=1e-3,call.rate=0.9)
## Imputation of the missing genotypes
data.imputed <- imputeSnpMatrix(data.scour$snpX, genes.info = data.scour$genes.info)
}
## Equivalent loading of the genotypes
load(system.file("extdata/dataImputed.Rdata", package="GeneGeneInteR"))

## Importation of the phenotype
resp <- system.file("extdata/response.txt", package="GeneGeneInteR")
Y  <- read.csv(resp, header=FALSE)

## estimation of the interaction between the 17 genes with the CLD method -- can take a few minutes
\dontrun{
GGI.res <- GGI(Y=Y, snpX=data.imputed$snpX, genes.info=data.imputed$genes.info,method="CLD")
}

## estimation of the interaction between 12 among the 17 genes with the default PCA method 
## Selection of 12 genes among 17
dta <- selectSnps(data.imputed$snpX, data.imputed$genes.info, c("bub3","CDSN","Gc","GLRX",
                  "PADI1","PADI2","PADI4","PADI6","PRKD3","PSORS1C1","SERPINA1","SORBS1"))
GGI.res <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info,method="PCA")

}
\references{
M. Emily. AGGrEGATOr: A Gene-based GEne-Gene interActTiOn test for case-control association studies, Statistical Application in Genetics and Molecular Biology, 15(2): 151-171, 2016.\cr
J. Li et al. Identification of gene-gene interaction using principal components. BMC Proceedings, 3 (Suppl. 7): S78, 2009.\cr
Qianqian Peng, Jinghua Zhao, and Fuzhong Xue. A gene-based method for detecting gene-gene co-association in a case-control study. European Journal of Human Genetics, 18(5) :582-587, 2010.\cr
Yuan, Z. et al. (2012): Detection for gene-gene co-association via kernel canonical correlation analysis, BMC Genetics, 13, 83.\cr
Larson, N. B. et al. (2013): A kernel regression approach to gene-gene interaction detection for case-control studies, Genetic Epidemiology, 37, 695-703.\cr
Indika Rajapakse, Michael D. Perlman, Paul J. Martin, John A. Hansen, and Charles Kooperberg. Multivariate detection of gene-gene interactions. Genetic Epidemiology, 36(6):622-630, 2012.\cr
X. Zhang et al. A PLSPM-based test statistic for detecting gene-gene co-association in genome-wide association study with case-control design. PLoS ONE, 8(4):e62129, 2013.\cr
J. Li, et al.. A gene-based information gain method for detecting gene-gene interactions in case-control studies. European Journal of Human Genetics, 23 :1566-1572, 2015.\cr
M.X. Li et al. GATES: A Rapid and Powerful Gene-Based Association Test Using Extended Simes Procedure, American Journal of Human Genetics, 88(3): 283-293, 2011.\cr
B. Jiang, X. Zhang, Y. Zuo and G. Kang. A powerful truncated tail strength method for testing multiple null hypotheses in one dataset. Journal of Theoretical Biology 277: 67-73, 2011.\cr
D.V. Zaykin, L.A. Zhivotovsky, P.H. Westfall and B.S. Weir. Truncated product method for combining P-values. Genetic epidemiology 22: 170-185, 2002.\cr
}


\seealso{
\code{\link{PCA.test}},
 \code{\link{CCA.test}}, \code{\link{KCCA.test}}, \code{\link{CLD.test}},
 \code{\link{PLSPM.test}}, \code{\link{GBIGM.test}}, \code{\link{plot.GGInetwork}},
 \code{\link{minP.test}}, \code{\link{gates.test}}, \code{\link{tTS.test}},
 \code{\link{tProd.test}}, \code{\link{imputeSnpMatrix}}
}

