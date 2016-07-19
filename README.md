# GeneGeneInteR

To install and load the package in R

```ruby
library(devtools)
install_github("MathieuEmily/GeneGeneInteR")
library(GeneGeneInteR)
```

 
Importation of genotypes with ImportFile function
Supported format are pedfile, PLINK, VCF (4.0) file, or genotypes imputed by IMPUTE2.

```ruby
#### Example of ped format with 17 genes
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")
dta <- importFile(file=ped, snps=info, pos=posi, pos.sep="\t")

## Importation of the phenotype
resp <- system.file("extdata/response.txt", package="GeneGeneInteR")
Y  <- read.csv(resp, header=FALSE)
```

Prior to the statistical analysis, dataset can be modified by applying filters to the SNPs (snpMatrixScour function) or by imputing missing genotypes (imputeSnpMatrix function). A subset of genes can also be selected with the select.snps function.


```ruby
## Filtering of the data: SNPs with MAF < 0.05 or p.value for HWE < 1e-3 are removed. No filtering is applied regarding missing data (NA.rate=1).
dta <- snpMatrixScour(snpX=dta$snpX,genes.info=dta$genes.info,min.maf=0.05,min.eq=1e-3,NA.rate=1)
## Imputation of the missing genotypes
dta <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)
## Selection of a subset of 12 genes
dta <- select.snps(dta$snpX, dta$genes.info, c("bub3","CDSN","Gc","GLRX","PADI1","PADI2","PADI4","PADI6","PRKD3","PSORS1C1","SERPINA1","SORBS1"))
```

Gene-based gene-gene interaction analysis can be performed by testing each pair of genes in the datatset (function GGI). 10 methods are implemented in the GeneGeneInteR package to test a pair of genes: 
- 6 Gene-Gene multidimensional methods
	- Principal Components Analysis - **PCA**
	- Canonical Correlation Analysis - **CCA**
	- Kernel Canonical Correlation Analysis - **KCCA**
	- Composite Linkage Disequilibrium - **CLD**
	- Partial Least Square Path Modeling - **PLSPM**
	- Gene-Based Information Gain Method - **GBIGM**
- 4 Gene-Gene interaction methods based on SNP-SNP interaction testing:
	- Minimum p-value test - **minP**
	- Gene Association Test using Extended Simes procedure - **GATES**
	- Truncated Tail Strength test - **tTS** 
	- Truncated p-value Product test - **tProd**
.

```ruby
## Testing for all pair of genes with the CLD method
GGI.res <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info,method="CLD")
```
Visualization of the results can be performed through either a matrix display (GGI.plot) or a network output (draw.network).

```ruby
## Plot of the results with default values
plot(GGI.res)
## Plot of the results with a threshold and an ordering of the genes.
plot(GGI.res,threshold=0.1,hclust.order=TRUE)

## Example of network with default threshold 0.05
plot(GGI.res,method="network")
## Example of network with threshold 0.01 where genes with no interaction are not plotted (plot.nointer=FALSE)
plot(GGI.res,threshold=0.1,plot.nointer=FALSE,method="network")
```


