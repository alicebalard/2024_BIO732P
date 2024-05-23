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
countdata <- read.table("../data/raw_counts_matrix.txt", header = T, row.names = 1)
coldata <- read.csv("../data/metadata.csv")

## we will make 2 groups: north and south
coldata$group = "North"
coldata$group[coldata$Latitude < 35] = "South"

### Keep only the distant populations
coldata = rbind(coldata[coldata$Latitude <35,],
                coldata[coldata$Latitude >40,])
countdata = countdata[colnames(countdata) %in% coldata$RNAseq_Run_NCBI]

## Make matrices with NAs
nbrElem = nrow(countdata)*ncol(countdata)

set.seed(1234)
## 5% values missing:

# Get indices of random values to replace
random_indices <- sample(1:nbrElem, nbrElem * 0.05)

# replace with NA
countdata5 = as.matrix(countdata)
countdata5[random_indices] <- NA
countdata5 = as.data.frame(countdata5)

## 10% values missing:

# Get indices of random values to replace
random_indices <- sample(1:nbrElem, nbrElem * 0.1)

# replace with NA
countdata10 = as.matrix(countdata)
countdata10[random_indices] <- NA
countdata10 = as.data.frame(countdata10)

## 30% values missing:

# Get indices of random values to replace
random_indices <- sample(1:nbrElem, nbrElem * 0.3)

# replace with NA
countdata30 = as.matrix(countdata)
countdata30[random_indices] <- NA
countdata30 = as.data.frame(countdata30)

write.table(countdata5, "../data/countdata5.txt")
write.table(countdata10, "../data/countdata10.txt")
write.table(countdata30, "../data/countdata30.txt")

coldata$group %>% table

#########################################


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

write.csv(results_df, file = "../data/results_DE-TRUTH.csv", row.names = F)
