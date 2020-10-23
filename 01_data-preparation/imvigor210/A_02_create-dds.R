
# Description -------------------------------------------------------------

# Coerce clinical data, counts, and genes into an SE


# Prepare Workspace -------------------------------------------------------

library(pointr)
library(readr)
library(DESeq2)

token <- sharepoint_token()


# Read in Data ------------------------------------------------------------

clin <- read_tsv("data/imvigor210/clin.tsv")

counts <- sharepoint_get("https://livejohnshopkins.sharepoint.com/sites/GBCIStorage/Shared%20Documents/Datasets/ImVigor210/IMvigor210counts.txt", token)
counts <- read.delim(counts)

genes <- sharepoint_get("https://livejohnshopkins.sharepoint.com/sites/GBCIStorage/Shared%20Documents/Datasets/ImVigor210/IMvigor210geneannot.txt", token)
genes <- read.delim(genes)

dds <- DESeqDataSetFromMatrix(as.matrix(counts), clin, design = ~1, rowData = genes)


# Write Data --------------------------------------------------------------

write_rds(dds, "./data/imvigor210/dds.Rds")

