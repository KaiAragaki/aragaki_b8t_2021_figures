# Creates a manifest of files to download


# Prepare Workspace -------------------------------------------------------

library(tidyverse)          # Quality of life
library(GenomicDataCommons) # Getting data
library(TCGAutils)          # For conversion of UUIDs to TCGA IDs
library(readxl)             # For reading in supplementary data


# Read in Data ------------------------------------------------------------

ffpe_cases <- read_tsv("./data/tcga-luad/ffpe-cases.tsv")


# Define Functions --------------------------------------------------------

getDuplicates <- function(x){
        x$shortID[duplicated(x$shortID)]
}


# Create Manifest ---------------------------------------------------------

manifest_count <- files() %>%
        GenomicDataCommons::filter(cases.project.project_id == 'TCGA-LUAD') %>% 
        GenomicDataCommons::filter(type == 'gene_expression') %>%
        GenomicDataCommons::filter(analysis.workflow_type == 'HTSeq - Counts') %>%
        manifest() %>%
        # Creates a 'path' feature that DESeq uses for its sampleTable argument
        mutate(path = paste(.$id, .$filename, sep = "/"))


# Get TCGA Barcodes from UUIDs --------------------------------------------

uuids <- manifest_count$id %>%
        UUIDtoBarcode(from_type = "file_id") %>%
        dplyr::rename(submitter_id = associated_entities.entity_submitter_id)


# Match Counts to Clinical ------------------------------------------------

matched <- manifest_count %>%
        left_join(uuids, by = c("id" = "file_id")) %>%
        mutate(shortID = substring(submitter_id, 1, 12))


# Filter Data -------------------------------------------------------------

table(str_extract(uuids$submitter_id, "(?<=TCGA-.{2}-.{4}-).{2}"))

# Select only tumors, not normal tissue
matched <- matched %>%
        dplyr::filter(str_detect(submitter_id, "^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-0"))

# Remove FFPE Cases
sum(matched$submitter_id %in% ffpe_cases$Barcode)
which(matched$submitter_id %in% ffpe_cases$Barcode)
matched <- matched[-which(matched$submitter_id %in% ffpe_cases$Barcode),]


# For duplicate samples, select RNA with later plate number
# as per GDAC Firehouse protocol: 
# http://gdac.broadinstitute.org/runs/stddata__2016_01_28/samples_report/LUAD_Replicate_Samples.html

# Filter and Remove Duplicate RNA Portions
matched <- matched %>% 
        mutate(rna_sample = str_sub(submitter_id, 1, 25),
               rna_portion = str_sub(submitter_id, 1, 20)) %>%
        arrange(desc(rna_sample)) %>%
        mutate(is_duplicated = duplicated(rna_portion)) %>%
        dplyr::filter(!is_duplicated) %>%
        dplyr::select(-is_duplicated, -rna_portion, -rna_sample)

sum(duplicated(matched$shortID))
matched <- matched[which(!duplicated(matched$shortID)),]


# Write Manifest ----------------------------------------------------------

write_tsv(matched, "./data/tcga-luad/01_manifest.tsv")
