
# Description -------------------------------------------------------------

# Creates a manifest of files to download


# Prepare Workspace -------------------------------------------------------

library(dplyr)              # Quality of life
library(readr)              # Writing CSV
library(GenomicDataCommons) # Getting data
library(TCGAutils)          # For conversion of UUIDs to TCGA IDs
library(readxl)             # For reading in supplementary data
library(jsonlite)


# Define Functions --------------------------------------------------------

getDuplicates <- function(x){
        x$shortID[duplicated(x$shortID)]
}


# Create Manifest ---------------------------------------------------------

geManifest <- files() %>%
        GenomicDataCommons::filter(cases.project.project_id == 'TCGA-HNSC') %>% 
        GenomicDataCommons::filter(type == 'gene_expression') %>%
        GenomicDataCommons::filter(analysis.workflow_type == 'HTSeq - Counts') %>%
        manifest() %>%
        # Creates a 'path' feature that DESeq uses for its sampleTable argument
        mutate(path = paste(.$id, .$filename, sep = "/"))


download.file("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5687509/bin/NIHMS911030-supplement-8.xlsx", 
              destfile = "./data/tcga-hnsc/sup.xlsx", method = "curl")

supData <- 
        read_xlsx("./data/tcga-hnsc/sup.xlsx", sheet = 2) %>%
        dplyr::filter(`RNA Seq` == "Yes") %>%
        dplyr::select(`Case ID`, `mRNA cluster`)

##
supData3 <- fromJSON("./data/tcga-hnsc/cases.2021-01-26.json")
download.file("http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/HNSC/20160128/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0.tar.gz", 
              destfile = "./data/tcga-hnsc/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0.tar.gz")

untar("./data/tcga-hnsc/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0.tar.gz", exdir = "./data/tcga-hnsc/")

clinical <- read_tsv("./data/tcga-hnsc/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0/HNSC.clin.merged.txt", col_names = F) %>% 
        t()
colnames(clinical) <- clinical[1,]

clinical <- clinical %>% 
        as_tibble()

clinical <- clinical[-1,]
supData2=clinical

untar("./data/tcga-hnsc/clinical.cases_selection.2021-01-26.tar.gz", exdir = "./data/tcga-hnsc/")
clinical2 <- read_tsv("./data/tcga-hnsc/clinical.tsv", col_names = T)
exposure <- read_tsv("./data/tcga-hnsc/exposure.tsv", col_names = T)
fhistory <- read_tsv("./data/tcga-hnsc/family_history.tsv", col_names = T)
##

uuids <- geManifest$id %>%
        UUIDtoBarcode(from_type = "file_id") %>%
        dplyr::rename(submitter_id = associated_entities.entity_submitter_id)

matched <- geManifest %>%
        left_join(uuids, by = c("id" = "file_id")) %>%
        # Creates a 'shortID' feature - this ID refers to the patient it came from
        # Various filtering methods will follow to ensure that there is only
        # one sample per patient - in its current state there are many duplicates
        mutate(shortID = substring(submitter_id, 1, 12)) %>%
        right_join(supData2, by = c("shortID" = "patient.bcr_patient_barcode"))


# Filter Data -------------------------------------------------------------

tcgaBarcode <- "^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-"
dupes <- getDuplicates(matched)

geManifest <- matched %>%
        # Selects for 'TCGA-XX-XXXX-01..." samples (vs 'TCGA-XX-XXXX-11...')
        # The '01' samples are tumor tissues, the '11' are normal tissues.
        dplyr::filter(grepl(paste0(tcgaBarcode, "01"),submitter_id))

dupes <- getDuplicates(geManifest)
geManifest <- geManifest %>%
        # Selects *against* duplicates that have the barcode "TCGA-XX-XXXX-01B"
        # These tissues were FFPE instead of frozen like everything else
        dplyr::filter(!(grepl(paste0(tcgaBarcode, "01B"), submitter_id) & (shortID %in% dupes)))

dupes <- getDuplicates(geManifest)
geManifest <- geManifest %>%
        # Removes one of each duplicate sample. These samples were removed in 
        # the original study.
        dplyr::filter(grepl(paste0(tcgaBarcode, "01A-11R-A277"), submitter_id) | !(shortID %in% dupes))


# Write Manifest ----------------------------------------------------------

write_csv(geManifest, "./data/tcga-hnsc/H_01_manifest.csv")