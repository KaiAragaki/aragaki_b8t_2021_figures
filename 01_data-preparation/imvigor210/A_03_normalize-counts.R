
# Description -------------------------------------------------------------

# Normalize raw counts and add as assay to dds object


# Prepare Workspace -------------------------------------------------------

library(DESeq2)
library(readr)
library(dplyr)


# Read in Data ------------------------------------------------------------

dds <- read_rds("./data/imvigor210/dds.Rds")
sigs <- read_rds("./data/signatures/signatures.Rds")

# Normalize ---------------------------------------------------------------

norm_counts <- dds %>%
        estimateSizeFactors() %>%
        vst() %>%
        assay()

assay(dds, 2) <- norm_counts
assayNames(dds)[[2]] <- "vst"

rownames(dds) <- make.unique(rowData(dds)$Symbol)

# Write Data --------------------------------------------------------------

write_rds(dds, "./data/imvigor210/norm_dds.Rds")
