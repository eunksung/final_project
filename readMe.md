TRGN510 Final Project
================

# **Title**

Differential Gene Expression in TCGA Cohort of Esophageal
Adenocarcinomas comparing by Alcohol Consumption with using DeSEQ2

# **Author**

Eun K. Sung

# **Overview of Project**

I am going to identify different gene expressions between esophageal
adenocarcionmas patients with alcohol consumption history and
non-alcohol consumption history. I utilize GDC Data Portal and have
found 86 samples fit within my cohort. I have picked 30 samples for
alcohol consumption cohort and 28 samples for non-alcohol consumption
cohort. I will utilize the package DeSEQ2 for this differential gene
expression analysis and follow the specific vignette: [Analyzing RNA-Seq
Data with
DeSEQ2](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

# **Data**

I will use the data from [GDC Data
Portal](https://portal.gdc.cancer.gov/repository) There are 1,138
samples of esophageal adenocarcinomas. By examining clinical data, there
are 86 samples defined by alcohol consumption history. All 28 samples of
non-alcohol consumption history have gene counts file by STAR. I have
selected 30 samples of alcohol consumption history with gene counts file
by STAR.

# **Method**
I will utilize unstrand counts in this differential expression analysis because I do not know of the strand-specific protocol of this dataset

Extract the column #4 which is unstranded counts of the sample
```{bash}
awk'{print $4}' AASX.tsv > AASX.txt
```

Remove the first line of the file without printing
```{bash}
tail -n +2 AASX.txt > AASX.tmp && mv AASX.tmp AASX.txt
```

Replace_header.sh to change the column of each .txt file to its file name
```{bash}
./replace_header.sh
```

Extract the gene_id, the column #1, from one of .tsv file
```{bash}
awk '{print $1}' A4OS.tsv > gene_id.txt
```

Remove the first line of the file without printing
```{bash}
tail -n +2 gene_id.txt > gene_id.tmp && mv gene_id.tmp gene_id.txt
```

Remove row #2 to #4
```{bash}
awk '!/^N_*/' gene_id.txt > gene_id.tmp && mv gene_id.tmp gene_id.txt
```

Merged counts
```{bash}
paste gene_id.txt A4OF.txt A4OJ.txt A4OR.txt A4OS.txt A4QS.txt A6BV.txt A6DN.txt A6DQ.txt A6FB.txt A6FH.txt A6FW.txt A6KZ.txt A6L4.txt A6L6.txt A6XG.txt A6Y0.txt A7BO.txt A7RE.txt A88T.txt A88V.txt A891.txt A8EQ.txt A8NF.txt A8NG.txt A8NH.txt A8NI.txt A8NJ.txt A8NL.txt A8NM.txt A8NR.txt A8NS.txt A8NU.txt A8NV.txt A8W8.txt A8WC.txt A8WG.txt A939.txt A93C.txt A93D.txt A93E.txt A9CJ.txt A9GF.txt A9GH.txt A9GI.txt A9GJ.txt A9GK.txt A9GL.txt A9GM.txt A9GN.txt A9GO.txt A9GQ.txt A9GR.txt A9W5.txt AA4D.txt AASW.txt AASX.txt > merged_gene_counts.txt
```

## Creating Sample Sheet

I change the sample id in the sample_sheet.tsv in order to match with the sample id in the count_matrix.txt. Also, I change the 1st column name to "sample_id" from "case_submitter_id"

```{bash}
sed 's/^TCGA-.*-//g' sample_sheet.tsv > sample_sheet.tmp && mv sample_sheet.tmp sample_sheet.tsv

sed 's/case_submitter_id/sample_id/g' sample_sheet.tsv > sample_sheet.tmp && mv sample_sheet.tmp sample_sheet.tsv
```

## Creating Sample Sheet

I change the sample id in the sample_sheet.tsv in order to match with the sample id in the count_matrix.txt. Also, I change the 1st column name to "sample_id" from "case_submitter_id"

---



## Load package into your library in RStudio

we are going to utilize these packages through this differential
expression analysis.

``` r
library(tibble)
library(tidyverse)
library(apeglm)
library(ggplot2)
library(vsn)
library(pheatmap)
library(ReportingTools)
library(DESeq2)
```

If you do not have any one of them, please use the command:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

    ## Bioconductor version '3.15' is out-of-date; the current release version '3.16'
    ##   is available with R version '4.2'; see https://bioconductor.org/install

``` r
library(BiocManager)
BiocManager::install("Named The Package You Need")
```

    ## Bioconductor version 3.15 (BiocManager 1.30.19), R 4.2.1 (2022-06-23)

    ## Installing package(s) 'Named The Package You Need'

    ## Warning: package 'Named The Package You Need' is not available for Bioconductor version '3.15'
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

    ## Old packages: 'bit', 'ggpubr', 'markdown', 'pbapply', 'pkgload', 'vctrs'

## Setting Working Directory and Wrangling the Raw Count Matrix and Sample Sheet

``` r
setwd('~/Desktop/final_project')

# Read in the matrix
count_matrix <- read.delim("~/Desktop/raw_data/merged_gene_counts.txt", header=T, sep="\t")

# Found Issue: I was not able to change gene_id (1st column) into a row name and found that I did not remove row #2 to #4 completely in the previous step
# I delete last 4 rows here
count_matrix <- head(count_matrix, - 4)

# Change column #1 into row name and delete the column #1 after that
row.names(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[-c(1)]

# Read in the sample sheet
sampletable <- read_tsv('~/Desktop/raw_data/sample_sheet.tsv')
```

``` r
# Change column #1 (sample_id) into row name
row.names(sampletable) <- sampletable$sample_id
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
# Change data type from character to factor
sampletable$alcohol_history <- as.factor(sampletable$alcohol_history)
```

## Create DESeq2 object

``` r
DES_dataset <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = sampletable,
                                         design = ~ alcohol_history)
```

## Filtering

``` r
# Number of gene before filtering
nrow(DES_dataset)
```

    ## [1] 60660

``` r
# Filtering to keep only rows that have at least 10 reads total
DES_dataset <- DES_dataset[rowSums(counts(DES_dataset)) > 10, ]

# Number of gene after filtering
nrow(DES_dataset)
```

    ## [1] 49404

# Performing standard differential expression analysis

``` r
DES_dataset <- DESeq(DES_dataset)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 3885 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

# Using “ReportingTools” to get result table of differential expression analysis

``` r
DES2Report <- HTMLReport(shortName = 'RNAseq_Analysis_with_DEseq2', title = 'Differential Expression Analysis in EAC', reportDirectory = "./reports")
publish(DES_dataset,DES2Report, pvalueCutoff=0.05, annotation.db="org.Mm.eg.db", factor = colData(DES_dataset)$alcohol_history, reportDir="./reports")
finish(DES2Report)
```

    ## [1] "./reports/RNAseq_Analysis_with_DEseq2.html"

This result table is more interactive compare to the regular result
table, which we will make it in the next step.  
You can double click “RNAseq_Analysis_with_DEseq2.html” file to open it.

# Generating a result table and print it

``` r
result_table <- results(DES_dataset)
result_table
```

    ## log2 fold change (MLE): alcohol history Yes vs No 
    ## Wald test p-value: alcohol history Yes vs No 
    ## DataFrame with 49404 rows and 6 columns
    ##                      baseMean log2FoldChange     lfcSE      stat     pvalue
    ##                     <numeric>      <numeric> <numeric> <numeric>  <numeric>
    ## ENSG00000000003.15 18010211.1      -0.552435 0.2672647  -2.06699 0.03873466
    ## ENSG00000000005.6  12377280.8      -0.778631 0.2427358  -3.20773 0.00133787
    ## ENSG00000000419.13 47198728.3      -0.153323 0.0493696  -3.10561 0.00189885
    ## ENSG00000000457.14  1732067.9      -0.177095 0.0746790  -2.37141 0.01772023
    ## ENSG00000000460.17     3241.2      -0.596481 0.2180447  -2.73559 0.00622683
    ## ...                       ...            ...       ...       ...        ...
    ## ENSG00000288662.1    19.43449      -0.515415  0.370249 -1.392078 0.16389891
    ## ENSG00000288665.1    96.27566      -0.933947  0.344594 -2.710283 0.00672258
    ## ENSG00000288669.1     3.71414      -0.284304  0.403720 -0.704211 0.48130143
    ## ENSG00000288670.1    59.23608       0.221601  0.153370  1.444883 0.14849083
    ## ENSG00000288674.1     2.60891       0.190368  0.352936  0.539384 0.58962197
    ##                         padj
    ##                    <numeric>
    ## ENSG00000000003.15 0.2480129
    ## ENSG00000000005.6  0.0363048
    ## ENSG00000000419.13 0.0451023
    ## ENSG00000000457.14 0.1650445
    ## ENSG00000000460.17 0.0923543
    ## ...                      ...
    ## ENSG00000288662.1  0.4962559
    ## ENSG00000288665.1  0.0967376
    ## ENSG00000288669.1  0.7817954
    ## ENSG00000288670.1  0.4724246
    ## ENSG00000288674.1  0.8420168

# Utilize “lfcShrink” function to shrink the effect size with apeglm method and print it

``` r
resultLFC <- lfcShrink(DES_dataset, coef = "alcohol_history_Yes_vs_No", type = "apeglm")
```

    ## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
    ##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
    ##     sequence count data: removing the noise and preserving large differences.
    ##     Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
resultLFC
```

    ## log2 fold change (MAP): alcohol history Yes vs No 
    ## Wald test p-value: alcohol history Yes vs No 
    ## DataFrame with 49404 rows and 5 columns
    ##                      baseMean log2FoldChange      lfcSE     pvalue      padj
    ##                     <numeric>      <numeric>  <numeric>  <numeric> <numeric>
    ## ENSG00000000003.15 18010211.1   -8.03881e-06 0.00144269 0.03873466 0.2480129
    ## ENSG00000000005.6  12377280.8   -5.71677e-01 0.30221611 0.00133787 0.0363048
    ## ENSG00000000419.13 47198728.3   -6.57038e-05 0.00144282 0.00189885 0.0451023
    ## ENSG00000000457.14  1732067.9   -3.31573e-05 0.00144262 0.01772023 0.1650445
    ## ENSG00000000460.17     3241.2   -1.30534e-05 0.00144269 0.00622683 0.0923543
    ## ...                       ...            ...        ...        ...       ...
    ## ENSG00000288662.1    19.43449   -3.91680e-06 0.00144269 0.16389891 0.4962559
    ## ENSG00000288665.1    96.27566   -8.18507e-06 0.00144269 0.00672258 0.0967376
    ## ENSG00000288669.1     3.71414   -1.80234e-06 0.00144269 0.48130143 0.7817954
    ## ENSG00000288670.1    59.23608    6.75892e-06 0.00144264 0.14849083 0.4724246
    ## ENSG00000288674.1     2.60891    2.54374e-06 0.00144268 0.58962197 0.8420168

# MA-plot to visualize the log2 fold change attributable to a given variable over the mean of normalized counts in the samples in the DES_dataset

``` r
#I will make one plot with the normal data and the other one with shrink the effect size
plotMA(result_table, ylim=c(-2,2))
```
![unnamed-chunk-10-1](https://user-images.githubusercontent.com/113143969/202277051-2dd51dbe-9523-4452-84b7-f8afe278be6e.png)

``` r
plotMA(resultLFC, ylim=c(-2,2))
```
![unnamed-chunk-10-2](https://user-images.githubusercontent.com/113143969/202277213-96049c79-189f-4c8d-9727-fdc662bf0406.png)

Set alpha to 0.05 for p-value and print the summary

``` r
result05 <- results(DES_dataset, alpha = 0.05)
summary(result05)
```

    ## 
    ## out of 49393 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 579, 1.2%
    ## LFC < 0 (down)     : 1146, 2.3%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 10545, 21%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
sum(result05$padj < 0.05, na.rm = TRUE)
```

    ## [1] 1725

# Use “plotCounts” function to make a plot for the read counts of single gene across the groups

``` r
plotCounts(DES_dataset, gene = which.min(result05$padj), intgroup = "alcohol_history")
```
![unnamed-chunk-12-1](https://user-images.githubusercontent.com/113143969/202277312-3557ca21-e15a-4720-b67f-7447a4665e23.png)

# Extracting transformed values by using variance stabilizing transformation (VST) and regularized log transformation

``` r
vsd <- vst(DES_dataset, blind = FALSE)
rld <- rlog(DES_dataset, blind = FALSE)
```

    ## rlog() may take a long time with 50 or more samples,
    ## vst() is a much faster transformation

``` r
head(assay(vsd), 3)
```

    ##                        A4OF     A4OJ     A4OR     A4OS     A4QS     A6BV
    ## ENSG00000000003.15 24.39564 24.87249 26.10379 22.33569 24.19637 23.26655
    ## ENSG00000000005.6  25.03573 25.67248 26.60153 22.84589 24.27888 23.40401
    ## ENSG00000000419.13 25.77906 25.82041 26.15114 25.25497 25.83509 25.53666
    ##                        A6DN     A6DQ     A6FB     A6FH     A6FW     A6KZ
    ## ENSG00000000003.15 23.18450 22.82781 23.21067 23.72965 23.16660 21.74093
    ## ENSG00000000005.6  23.05322 23.12316 22.72474 22.49905 23.37450 22.70998
    ## ENSG00000000419.13 25.41788 25.34583 25.46201 25.33988 25.59348 25.20305
    ##                        A6L4     A6L6     A6XG     A6Y0     A7BO     A7RE
    ## ENSG00000000003.15 23.13067 26.86190 23.85458 24.26617 24.90104 24.23776
    ## ENSG00000000005.6  23.31683 23.72701 23.29407 23.69903 23.00900 22.77393
    ## ENSG00000000419.13 25.56698 25.85192 25.50339 25.75041 25.60723 25.52708
    ##                        A88T     A88V     A891     A8EQ     A8NF     A8NG
    ## ENSG00000000003.15 24.43949 24.46294 24.82121 24.54520 23.33364 23.15901
    ## ENSG00000000005.6  23.92320 23.26595 23.19533 23.23656 22.62250 22.56560
    ## ENSG00000000419.13 25.73173 25.52122 25.50800 25.55180 25.59322 25.26414
    ##                        A8NH     A8NI     A8NJ     A8NL     A8NM     A8NR
    ## ENSG00000000003.15 22.70305 22.75252 24.95803 22.75865 22.69078 23.48815
    ## ENSG00000000005.6  22.90518 22.82301 22.92601 22.57963 22.70722 23.04360
    ## ENSG00000000419.13 25.47422 25.24726 25.38143 25.35165 25.37954 25.40272
    ##                        A8NS     A8NU     A8NV     A8W8     A8WC     A8WG
    ## ENSG00000000003.15 23.72539 25.17372 24.61568 22.29872 23.29910 22.82587
    ## ENSG00000000005.6  22.51950 23.83963 22.73376 22.60873 22.65958 22.53612
    ## ENSG00000000419.13 25.33734 25.54606 25.36008 25.44332 25.29114 25.38539
    ##                        A939     A93C     A93D     A93E     A9CJ     A9GF
    ## ENSG00000000003.15 25.30743 22.60844 23.25500 24.20573 23.09432 23.85000
    ## ENSG00000000005.6  25.15437 23.05151 22.91225 22.80760 22.82416 22.91301
    ## ENSG00000000419.13 25.63279 25.37991 25.37526 25.25188 25.40529 25.32513
    ##                        A9GH     A9GI     A9GJ     A9GK     A9GL     A9GM
    ## ENSG00000000003.15 23.18483 24.39200 22.89887 24.31939 23.49959 22.97108
    ## ENSG00000000005.6  22.81907 23.60607 23.05208 23.43474 22.81993 22.87908
    ## ENSG00000000419.13 25.26924 25.53077 25.52100 25.57128 25.40068 25.40558
    ##                        A9GN     A9GO     A9GQ     A9GR     A9W5     AA4D
    ## ENSG00000000003.15 22.93516 23.56825 23.13055 24.22145 23.34592 24.23779
    ## ENSG00000000005.6  22.96641 23.40530 22.50144 23.06731 23.05800 22.99888
    ## ENSG00000000419.13 25.43851 25.53534 25.28868 25.39069 25.30318 25.44204
    ##                        AASW     AASX
    ## ENSG00000000003.15 23.71378 23.30935
    ## ENSG00000000005.6  23.24050 23.08809
    ## ENSG00000000419.13 25.54868 25.57369

# Effects of transformations on the variance

``` r
# this gives log2(n+1)
ntd <- normTransform(DES_dataset)
meanSdPlot(assay(ntd))
```

    ## Found more than one class "simpleUnit" in cache; using the first, from namespace 'ggbio'

    ## Also defined by 'hexbin'

![unnamed-chunk-14-1](https://user-images.githubusercontent.com/113143969/202277475-9ef93faa-067f-4d1a-9a95-2f558249c6fc.png)

``` r
# using VST
meanSdPlot(assay(vsd))
```

![unnamed-chunk-14-2](https://user-images.githubusercontent.com/113143969/202277529-200013e9-12bd-4143-a623-7f34266228a1.png)

``` r
# using rlog
meanSdPlot(assay(rld))
```

![unnamed-chunk-14-3](https://user-images.githubusercontent.com/113143969/202277597-1ff33078-4c7b-46e7-a8cb-b0ec6fd3908c.png)

# Data quality assessment by sample clustering and visualization

``` r
select <- order(rowMeans(counts(DES_dataset,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(DES_dataset)[,c("sample_id", "alcohol_history")])
row.names(df) <- df$sample_id
df <- df[-c(1)]
```

# Using normal transform, variance stabilizing transformation, and regularized log transformation to make a heatmap

``` r
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
```

    ## Found more than one class "unit" in cache; using the first, from namespace 'ggbio'

    ## Also defined by 'hexbin'

    ## Found more than one class "simpleUnit" in cache; using the first, from namespace 'ggbio'

    ## Also defined by 'hexbin'

![unnamed-chunk-16-1](https://user-images.githubusercontent.com/113143969/202277724-246dbff9-2b4e-468a-b04f-d236f7d8d8e7.png)

``` r
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
```

![unnamed-chunk-16-2](https://user-images.githubusercontent.com/113143969/202277766-06ff22a2-9ac7-41c6-a5c6-32916278dd7e.png)

``` r
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
```

![unnamed-chunk-16-3](https://user-images.githubusercontent.com/113143969/202277889-d0749f9c-b00d-4a48-ac53-5aeff1b07653.png)

#Sample-to-Sample distances

``` r
sampleDists <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(sampleDists)
pheatmap(DistMatrix)
```

![unnamed-chunk-17-1](https://user-images.githubusercontent.com/113143969/202277961-02830854-56bf-476c-b187-660cee9c610c.png)

# Principal Component Analysis Plot

``` r
plotPCA(vsd, intgroup="alcohol_history")
```

![unnamed-chunk-18-1](https://user-images.githubusercontent.com/113143969/202278003-3d889ce5-3e16-41cf-a63a-3ee6299d0c6a.png)
