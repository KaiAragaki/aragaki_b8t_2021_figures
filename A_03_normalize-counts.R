
# Description -------------------------------------------------------------

# This file will take the outputs of A_01_create-manifest.R and 
# B_01_tidy-clinical.R, merges them, adds HGCN identifiers, and normalizes counts.
# These files will be written as an .Rds
#BiocManager::install("DESeq2")


# Initialize Workspace ----------------------------------------------------

library(tidyverse)          # 
library(DESeq2)             # DEGS, clustering
library(readxl)             # For reading in normalized counts for HGCN ID
library(biomaRt)            # Get HGCN, ensembl, etc IDs
library(GenomicDataCommons) # Get data


# Download Files ----------------------------------------------------------

geManifest <- read_csv('H:/Rusty_B8T/data/TCGA_BLCA/A_01_manifest.csv') %>%
  mutate(shortID = tolower(shortID))

gdc_set_cache("H:/Rusty_B8T/data/gdcdata")
fnames <- lapply(geManifest$id, gdcdata)

View(fnames)
# Attach Tumor Calls and Clinical Data ------------------------------------

clin <- read_csv("H:/Rusty_B8T/data/TCGA_BLCA/B_01_tidy-clinical.csv")

tumorData <- geManifest %>%
  dplyr::select(c(1, 6, 8, 9)) %>% 
  left_join(clin, by = c("shortID" = "patient.bcr_patient_barcode")) %>% 
  mutate(shortID = toupper(shortID)) 
View(tumorData)

# Generate Sample Table ---------------------------------------------------

sampleTable <- data.frame(sampleName = tumorData$id,
                          fileName = tumorData$path)
View(sampleTable)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory   = gdc_cache(),
                                  design      = ~ 1)


# Strip Ensembl Version Number --------------------------------------------

rownames(dds) <- str_replace(rownames(dds), ".[0-9]+$", "")

View(dds)
# Get HGNC Symbol ---------------------------------------------------------

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ids <- getBM(attributes = c('hgnc_symbol','entrezgene_id', 'ensembl_gene_id'),
             filters = 'ensembl_gene_id',
             values = rownames(dds),
             mart = mart)

# Many (~5%, or 3000) transcripts mapped to now deprecated genes and thus return NA
ids2 <- ids[match(rownames(dds), ids$ensembl_gene_id),]
mcols(dds) <- cbind(mcols(dds), ids2)

samples <- tumorData[match(colnames(dds), tumorData$id),]
colData(dds) <- cbind(colData(dds), samples)
colnames(dds) <- colData(dds)$shortID
rownames(dds) <- rowData(dds)$hgnc_symbol


# Remove NA HUGO symbols --------------------------------------------------

dds <- dds[which(!is.na(rownames(dds))),]


# Normalize Counts use vst() instead of normtransform()log2_norm_counts--------------------------------------------------------

normCounts <- dds %>%
  estimateSizeFactors() %>%
  vst() %>%
  assay()



assay(dds, 2) <- normCounts
assayNames(dds)[[2]] <- "vst"

View(assay(dds,2))
norm_counts_tibble <- normCounts %>%
  as_tibble(rownames = "Genes")

View(norm_counts_tibble)
# Write Data --------------------------------------------------------------

write_rds(dds, "H:/Rusty_B8T/data/TCGA_BLCA/A_02_normalized-counts.Rds")
