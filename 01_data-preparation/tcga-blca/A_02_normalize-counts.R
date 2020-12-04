
# Description -------------------------------------------------------------

# Hi Jerry

# This file will take the outputs of A_01_create-manifest.R and 
# B_01_tidy-clinical.R, merges them, adds HGCN identifiers, and normalizes counts.
# These files will be written as an .Rds


# Initialize Workspace ----------------------------------------------------

library(tidyverse)          # Quality of Life
library(DESeq2)             # DEGS, clustering
library(readxl)             # For reading in normalized counts for HGCN ID
library(biomaRt)            # Get HGCN, ensembl, etc IDs
library(GenomicDataCommons) # Get data


# Download Files ----------------------------------------------------------

geManifest <- read_csv('./data/tcga-blca/A_01_manifest.csv') %>%
        mutate(shortID = tolower(shortID))
gdc_set_cache("./data/gdcdata")
fnames <- lapply(geManifest$id, gdcdata)


# Attach Tumor Calls and Clinical Data ------------------------------------

clin <- read_csv("./data/tcga-blca/B_01_tidy-clinical.csv")

tumorData <- geManifest %>%
        dplyr::select(c(1, 6, 8, 9)) %>% 
        left_join(clin, by = c("shortID" = "patient.bcr_patient_barcode")) %>% 
        mutate(shortID = toupper(shortID)) 


# Generate Sample Table ---------------------------------------------------

sampleTable <- data.frame(sampleName = tumorData$id,
                          fileName = tumorData$path)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory   = gdc_cache(),
                                  design      = ~ 1)


# Strip Ensembl Version Number --------------------------------------------

rownames(dds) <- str_replace(rownames(dds), ".[0-9]+$", "")


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


# Normalize Counts --------------------------------------------------------

normCounts <- dds %>%
        estimateSizeFactors() %>%
        normTransform() %>%
        assay()

assay(dds, 2) <- normCounts
assayNames(dds)[[2]] <- "log2_norm_counts"

norm_counts_tibble <- normCounts %>%
        as_tibble(rownames = "Genes")


# Write Data --------------------------------------------------------------

write_rds(dds, "./data/tcga-blca/A_02_normalized-counts.Rds")
