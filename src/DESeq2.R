# Differential gene expression analysis using DESeq2 package
# getwd()

# Load libraries

library(DESeq2)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(airway)

library(tidyverse)
library(pheatmap)

# Step 1: Preparing count data ----------------

## Import read counts data

counts_data = read.table("./data/airway.tsv")
head(counts_data)
str(counts_data)

## Import sample information

colData = read.table("./data/airway_colData.tsv")
head(colData)
str(colData)

## Make sure the row names in colData (sample names) matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# Are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: Construct a DESeq data set object ----------

dds = DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dex)

dds

## Data filtering: removing genes with low counts across all the samples
## keeping only those genes that have at least 10 reads total

keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

dds

## set the factor level
dds$dex = relevel(dds$dex, ref = "untrt")

## NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------

dds = DESeq(dds)

# Step 4: rlog transformation ---------------------------------------------

rld = rlog( dds )
assay(rld)[ 1:3, 1:3]
assay(dds)[ 1:3, 1:3]


# Step 5: QC based on normalized counts -------------------------------------------

## Sample distances

sampleDists = dist( t( assay(rld) ) )
as.matrix( sampleDists )[ 1:3, 1:3 ]

sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = paste( rld$dex, 
                                     rld$Run, sep="-" )
colnames(sampleDistMatrix) = NULL 
pheatmap(sampleDistMatrix, cutree_rows = 2, cutree_cols = 2)


## PCA Plot

plotPCA( rld, intgroup = c("dex"))

# Step 6: Statistical testing and exploratory analysis ----------------

res = results(dds)
res

summary(res)

res0.01 = results(dds, alpha = 0.01)
summary(res0.01)

## Assign gene symbols to resulted genes
?mapIds
columns(org.Hs.eg.db)
symbols = mapIds(org.Hs.eg.db, keys = rownames(res0.01), column = "SYMBOL", 
                  keytype = "ENSEMBL", multiVals = "first")

str(res0.01)

str(symbols)
identical(names(symbols), rownames(res0.01))
res0.01$gene = symbols

df_res0.01 = as.data.frame(res0.01)
sigDf = df_res0.01 %>% 
  filter(!is.na(gene), padj < 0.01, abs(log2FoldChange) > 2) %>% 
  arrange(desc(log2FoldChange))

# contrasts
resultsNames(dds)

# e.g.: treated_4hrs, treated_8hrs, untreated

# results(dds, contrast = c("dex", "treated_4hrs", "untrt"))
?results

# Step 7: Data visualization ------------------------------------------------------

# MA plot

plotMA(res)
maData = plotMA(res, returnData = TRUE)
plotMA( res, ylim = c(-4, 4), alpha = 0.01, colSig = "red" )

ggplot(data = maData) +
  aes(x = log2(mean), y = lfc, colour = isDE) +
  geom_point()
  
plotDispEsts( dds, ylim = c(1e-6, 1e1) )

hist( res$pvalue, breaks=20, col="grey" )
hist( res$padj, breaks=20, col="grey" )

## Volcano plot

EnhancedVolcano(res0.01,
                lab = res0.01$gene,
                x = 'log2FoldChange',
                y = 'padj')

## Heatmap

mat =  assay(rld)[ rownames(sigDf), ]
boxplot(mat)
mat[mat>14] = 14
boxplot(mat)
dim(mat)
pMat = head(mat, 60)
pheatmap(pMat, 
         treeheight_row = 0, 
         treeheight_col = 0, 
         gaps_row = c(4), 
         cutree_cols = 2, 
         labels_row = sigDf$gene, 
         annotation_col = colData %>% dplyr::select(dex),
         fontsize_row = 10,
         fontsize_col = 10)

# Step 8: Export list of up-regulated and down regulated genes --------------

str(sigDf)

## Up-regulated genes 

sigDf %>% filter(log2FoldChange > 2) %>% 
  select(gene) %>% 
  write.table(., file = "upGene.tsv", row.names = FALSE)

## Down-regulated genes 

sigDf %>% filter(log2FoldChange < 2) %>% 
  select(gene) %>% 
  write.table(., file = "downGene.tsv", row.names = FALSE)

