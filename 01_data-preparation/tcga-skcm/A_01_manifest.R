
# Description -------------------------------------------------------------

# Creates a manifest of files to download


# Prepare Workspace -------------------------------------------------------

library(tidyverse)          # Quality of life
library(GenomicDataCommons) # Getting data
library(TCGAutils)          # For conversion of UUIDs to TCGA IDs


# Create Manifest ---------------------------------------------------------

manifest_count <- files() %>%
        GenomicDataCommons::filter(cases.project.project_id == 'TCGA-SKCM'
                                   & type == 'gene_expression'
                                   & analysis.workflow_type == 'HTSeq - Counts') %>%
        manifest() %>%
        # Creates a 'path' feature that DESeq uses for its sampleTable argument
        mutate(path = paste(.$id, .$filename, sep = "/"))


# Get TCGA Barcodes from UUIDs --------------------------------------------

uuids <- manifest_count$id %>%
        UUIDtoBarcode(from_type = "file_id") %>%
        rename(submitter_id = associated_entities.entity_submitter_id)


# Match Counts to UUIDs ---------------------------------------------------

matched <- manifest_count %>%
        left_join(uuids, by = c("id" = "file_id")) %>%
        mutate(shortID = str_sub(submitter_id, 1, 12))


# Filter Data -------------------------------------------------------------

# Select only tumors, not normal tissue
matched <- matched %>%
        dplyr::filter(str_detect(submitter_id, "^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-0"))

# Look at types of samples
table(str_extract(matched$submitter_id, "(?<=^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-)[:digit:]{2}"))

# 01 = primary 
# 06 = metastasis
# 07 = additional metastasis

# There are roughly three patients that have multiple tumors to them. We should keep these in mind
# We should also keep in mind that the tumors we have are primary and metastases


# Write Manifest ----------------------------------------------------------

write_csv(matched, "./data/tcga-skcm/A_01_manifest.csv")
