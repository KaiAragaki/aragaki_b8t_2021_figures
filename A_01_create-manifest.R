
# Description -------------------------------------------------------------

# Creates a manifest of files to download

# supData: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5687509/bin/NIHMS911030-supplement-8.xlsx


# Prepare Workspace -------------------------------------------------------

library(dplyr)              # Quality of life
library(readr)              # Writing CSV
library(GenomicDataCommons) # Getting data
library(TCGAutils)          # For conversion of UUIDs to TCGA IDs
library(readxl)             # For reading in supplementary data


# Define Functions --------------------------------------------------------

getDuplicates <- function(x){
  x$shortID[duplicated(x$shortID)]
}


# Create Manifest ---------------------------------------------------------

geManifest <- files() %>%
  GenomicDataCommons::filter(cases.project.project_id == 'TCGA-BLCA') %>% 
  GenomicDataCommons::filter(type == 'gene_expression') %>%
  GenomicDataCommons::filter(analysis.workflow_type == 'HTSeq - Counts') %>%
  manifest() %>%
  # Creates a 'path' feature that DESeq uses for its sampleTable argument
  mutate(path = paste(.$id, .$filename, sep = "/"))

##
#write_csv(geManifest, "H:/Rusty_B8T/data/TCGA_BLCA/manifest_tcga.csv")
##
#View(geManifest)

supData <- 
  read_xlsx("H:/Rusty_B8T/data/TCGA_BLCA/NIHMS911030-supplement-8.xlsx", 
            sheet = 2) %>%
  dplyr::filter(`RNA Seq` == "Yes") %>%
  dplyr::select(`Case ID`, `mRNA cluster`)

#View(supData)

uuids <- geManifest$id %>%
  UUIDtoBarcode(from_type = "file_id") %>%
  dplyr::rename(submitter_id = associated_entities.entity_submitter_id)

#View(uuids)

matched <- geManifest %>%
  left_join(uuids, by = c("id" = "file_id")) %>%
  # Creates a 'shortID' feature - this ID refers to the patient it came from
  # Various filtering methods will follow to ensure that there is only
  # one sample per patient - in its current state there are many duplicates
  mutate(shortID = substring(submitter_id, 1, 12)) %>%
  right_join(supData, by = c("shortID" = "Case ID"))

#View(matched)

# Filter Data -------------------------------------------------------------

tcgaBarcode <- "^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-"


dupes <- getDuplicates(matched)
View(sort(dupes))

geManifest <- matched %>%
  # Selects for 'TCGA-XX-XXXX-01..." samples (vs 'TCGA-XX-XXXX-11...')
  # The '01' samples are tumor tissues, the '11' are normal tissues.
  dplyr::filter(grepl(paste0(tcgaBarcode, "01"),submitter_id))

dupes <- getDuplicates(geManifest)


geManifest <- geManifest %>%
  # Selects *against* duplicates that have the barcode "TCGA-XX-XXXX-01B"
  # These tissues were FFPE instead of frozen like everything else
  dplyr::filter(!(grepl(paste0(tcgaBarcode, "01B"), submitter_id) & (shortID %in% dupes)))

View(geManifest)


dupes <- getDuplicates(geManifest)
View(dupes)
geManifest <- geManifest %>%
  # Removes one of each duplicate sample. These samples were removed in 
  # the original study.? what is the original study?
  dplyr::filter(grepl(paste0(tcgaBarcode, "01A-11R-A277"), submitter_id) | !(shortID %in% dupes))

View(geManifest)

# Write Manifest ----------------------------------------------------------

write_csv(geManifest, "H:/Rusty_B8T/data/TCGA_BLCA/A_01_manifest.csv")