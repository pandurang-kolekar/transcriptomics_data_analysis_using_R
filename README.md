# Transcriptomics data analysis using R/Bioconductor

## Abstract

RNA sequencing (RNASeq) has become a method of choice for transcriptome profiling, however the analysis of the massive amount of data generated by large‐scale RNASq still remains a challenge. Typically RNA‐seq data analyses consist of (1) alignment of millions of short sequencing reads to a reference genome or *de novo* assembled transcriptome, including the identification of splicing events; (2) quantification of expression levels of genes, transcripts, and exons; (3) differential analysis of gene expression among different biological conditions; and (4) biological interpretation of those differentially expressed genes.

Various Bioinformatics pipeline exist for quantification of expression levels of genes/transcripts, which ultimately generate a matrix of read counts for genes/transcripts in given samples. This count matrix is then subsequently used for differential expression analysis using R/Bioconductor packages.

The workshop session will demonstrate the applications of R/Bioconductor packages for differential expression analysis and visualization of transcriptomics data. During this workshop, participants will learn how to import count data, pre-process and normalize data, perform statistical analyses to find differentially expressed genes and generate publication ready figures to report the results.

## Software Installation

The transcriptomics data analysis will be carried out using [R (version 4.2.1)](https://cloud.r-project.org/) and [Bioconductor (version 3.15)](https://bioconductor.org/install/) packages.

We will be using [RStudio](https://www.rstudio.com/products/rstudio/), an integrated development environment (IDE) for R to edit & execute scripts and visualize plots.

First install [R (version 4.2.1)](https://cloud.r-project.org/) and then [RStudio](https://www.rstudio.com/products/rstudio/). Once the installation is completed, run the RStudio.

***Overview of RStudio IDE***

[Youtube Video](https://youtu.be/n3uue28FD0w)

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/n3uue28FD0w/0.jpg)](https://www.youtube.com/watch?v=n3uue28FD0w)

### Install R packages

Run RStudio and execute following commands on R Console to install R packages. [Tidyverse](https://www.tidyverse.org/) is a collection of standard R packages that are widely used in data transformation and visualization.

```{r eval=FALSE, include=TRUE}
install.packages(c("tidyverse", "pheatmap"), dependencies = TRUE)
```

### Install Bioconductor

Run RStudio and execute following commands on R Console to install the [Bioconductor](https://bioconductor.org/install/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and related packages.

```{r eval=FALSE, include=TRUE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")
BiocManager::install(c("DESeq2", "edgeR", "org.Hs.eg.db", "EnhancedVolcano"))
```

### Test installations

```{r eval=FALSE, include=TRUE}
library(tidyverse)
library(pheatmap)
library(DESeq2)
library(edgeR)
library(org.Hs.eg.db)
library(EnhancedVolcano)
```

## Data Set

Bioinformatics pipelines designed for transcriptomics or RNASeq data analysis generate a matrix of read counts for genes/transcripts in given samples. The rows in matrix represent either gene or transcript IDs and columns represent samples. The values in the matrix correspond to number of reads aligned to corresponding gene or transcript in related sample.

We will be using previously published "***airway***" data set, which was described in the following publication,

Himes BE, Jiang X, Wagner P, Hu R, Wang Q, et al. (2014) *RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells*. ***PLOS ONE*** 9(6): e99625. <https://doi.org/10.1371/journal.pone.0099625>

**PMID**: [24926665](https://pubmed.ncbi.nlm.nih.gov/24926665/). **GEO**: [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).

Authors used RNASeq experiment to characterize the human airway smooth muscle (HASM) transcriptome at baseline and under dexamethasone (DEX) asthma treatment. RNASeq data from HASM cell lines with untreated (n=4) and treated (n=4) samples were processed through Bioinformatics pipeline as described in the paper to generate read count matrix across genes and samples.

-   The count matrix data set can be downloaded either in the form of Bioconductor [*airway*](https://bioconductor.org/packages/release/data/experiment/vignettes/airway/inst/doc/airway.html) package.

```{r eval=FALSE, include=TRUE}
BiocManager::install("airway")
library("airway")
data(airway)
dim(airway)
str(assay(airway))
colData(airway)
```

-   Tab-delimited text files for data sets are saved in the *"./data"* directory of this repository

    ```{r}
    list.files("./data")
    ```

-   Data set can also be downloaded as tab-delimited files from [here](https://drive.google.com/drive/folders/1eosGdWh0zzEXLg1v61BpMQAk60YpoYJ7?usp=sharing) (Download both the files from Google drive, [airway.tsv](https://drive.google.com/file/d/11dujOxITpI3KHujyBTojfNBueLXbAIzk/view?usp=sharing) and [airway_colData.tsv](https://drive.google.com/file/d/1Vji3F2M4VUIoQOV_4_yeRGdPNhIjU1sJ/view?usp=sharing)).

## Source Code for Hands-on session

An R script with commands for step wise differential expression analysis using DESeq2 package is saved in the *"./src"* directory.

```{r}
"./src/DESeq2.R"
```
