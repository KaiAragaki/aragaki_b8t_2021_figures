
# Description -------------------------------------------------------------

# Prepare Riaz melanoma data
# Downloaded from here: 
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE91061&format=file&file=GSE91061%5FBMS038109Sample%2Ehg19KnownGene%2Eraw%2Ecsv%2Egz


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(DESeq2)
library(biomaRt)


# Read in Data ------------------------------------------------------------

counts <- read_csv("./data/riaz/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv")
clinical <- read_rds("./data/riaz/tidy_clinical.Rds")


# Prepare rowData ---------------------------------------------------------

row_data <- DataFrame(entrez = counts$X1)


# Prepare colData ---------------------------------------------------------

col_data <- tibble(id = colnames(counts)[-1]) %>% 
        left_join(clinical, by = "id")


# Prepare SummarizedExperiment --------------------------------------------

se <- SummarizedExperiment(as.matrix(counts[-1]), rowData = row_data, colData = col_data)
dds <- DESeqDataSet(se, design = ~1)
        
norm_counts <- dds %>%
        estimateSizeFactors() %>%
        vst() %>%
        assay()

assay(dds, 2) <- norm_counts
assayNames(dds)[[2]] <- "vst"

norm_counts_tibble <- norm_counts %>%
        as_tibble(rownames = "genes")

rownames(dds) <- rowData(dds)$entrez

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ids <- getBM(attributes = c('hgnc_symbol','entrezgene_id', 'ensembl_gene_id'),
             filters = 'entrezgene_id',
             values = rownames(dds),
             mart = mart)

joined <- left_join(as_tibble(rowData(dds)), as_tibble(ids), by = c("entrez" = "entrezgene_id"))

joined <- joined[isUnique(joined$entrez), ]

joined <- joined[!is.na(joined$hgnc_symbol) & joined$hgnc_symbol != "", ]

dds <- dds[which(rownames(dds) %in% joined$entrez), ]

rowData(dds) <-  DataFrame(left_join(as_tibble(rowData(dds)), joined, by = "entrez"))

rownames(dds) <- rowData(dds)$hgnc_symbol

# Write Summarized Experiment ---------------------------------------------

write_rds(dds, "./data/riaz/dds.Rds", compress = "gz")
