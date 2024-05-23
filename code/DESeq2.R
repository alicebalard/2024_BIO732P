## Differential expression analysis based on the data from Mack et al. 2018
## A. Balard, 2024

## Source: https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html

## Install BiocManager if not yet
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

library('DESeq2') # Package for DE analysis
library(ggplot2) # To plot data
library(dplyr)
library(tibble)
library(gprofiler2)
## BiocManager::install('EnhancedVolcano') # package for pretty volcano plots
library(EnhancedVolcano)

## For the differential expression analysis, we need:
# countdata: a table with the fragment counts
# coldata: a table with information about the samples
countdata <- read.table("../data/countdata5.txt", header = T, row.names = 1)
coldata <- read.csv("../data/metadata.csv")

## we will make 2 groups: north and south
coldata$group = "North"
coldata$group[coldata$Latitude < 35] = "South"

### Keep only the distant populations
coldata = rbind(coldata[coldata$Latitude <35,],
                coldata[coldata$Latitude >40,])
countdata = countdata[colnames(countdata) %in% coldata$RNAseq_Run_NCBI]

# NB: it is very important to check manually that the columns of the count matrix
# correspond to the rows of the sample information table.
colnames(countdata) == coldata$RNAseq_Run_NCBI
coldata <- coldata[match(colnames(countdata), coldata$RNAseq_Run_NCBI), ]
colnames(countdata) == coldata$RNAseq_Run_NCBI # now we can proceed

## rm NA -> complete cases
countdata = na.omit(countdata)

## Create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ Sex + group) # we want to test for group while controlling for covariate sex

## Run analysis (everything from normalization to linear modeling is included in this function)
dds <- DESeq(dds)
res <- results(dds)

## Summarize results
summary(res)

### Extracting significant differentially expressed genes

## Convert Ensembl names to gene names (if they are known)
results_df = as.data.frame(res)
results_df$Ensembl_id = row.names(results_df)
results_df = results_df[order(results_df$padj),]

results_genes = gconvert(row.names(res), organism = "mmusculus",
                         target = "ENTREZGENE_ACC", filter_na = FALSE)
# add the gene names
results_df = merge(results_df,
                   results_genes[,c("input", "target", "name", "description")],
                   by.x = "Ensembl_id", by.y = "input")

results_df$Name <- ifelse(is.na(results_df$name), results_df$Ensembl_id, results_df$name)

# Subset the results to keep only significant genes
results_df[results_df$padj < 0.05 & !is.na(results_df$padj),]

## Volcano plot
EnhancedVolcano(results_df,
    lab = results_df$Name,
    x = 'log2FoldChange',
    y = 'padj',
    drawConnectors = TRUE)

## Final result
# write.csv(results_df, file = "../data/results_DE-5%.csv", row.names = F)
