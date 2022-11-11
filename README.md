<h1>Title</h1>
Differential Gene Expression in TCGA Cohort of Esophageal Adenocarcinomas comparing by Alcohol Consumption with using DeSEQ2

<h1>Author</h1>
Eun K. Sung

<h1>Overview of Project</h1>
I am going to identify different gene expressions between esophageal adenocarcionmas patients with alcohol consumption history and non-alcohol consumption history. I utilize GDC Data Portal and have found 86 samples fit within my cohort. I have picked 30 samples for alcohol consumption cohort and 28 samples for non-alcohol consumption cohort. I will utilize the package DeSEQ2 for this differential gene expression analysis and follow the specific vignette: <a href="http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html">Analyzing RNA-Seq Data with DeSEQ2</a>

<h1>Data</h1>
I will use the data from <a href="https://portal.gdc.cancer.gov/repository">GDC Data Portal</a>. There are 1,138 samples of esophageal adenocarcinomas. By examining clinical data, there are 86 samples defined by alcohol consumption history.  All 28 samples of non-alcohol consumption history have gene counts file by STAR. I have selected 30 samples of alcohol consumption history with gene counts file by STAR.

