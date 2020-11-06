
# Description -------------------------------------------------------------

# Prepare van Allen melanoma data
# Clinical data downloaded from here: https://www.cbioportal.org/study/summary?id=skcm_dfci_2015


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(SummarizedExperiment)
library(biomaRt)


# Read in Data ------------------------------------------------------------

rna <- read_tsv("./data/van-allen/data_RNA_Seq_expression_median.txt") %>% 
        column_to_rownames('Entrez_Gene_Id') %>% 
        as.matrix()
        
clin <- read_tsv("./data/van-allen/data_clinical_patient.txt", col_names = F, skip = 5)
clin_colnames <- read_tsv("./data/van-allen/data_clinical_patient.txt", n_max = 1, col_names = F)
clin_colnames <- tolower(clin_colnames) %>% 
        str_replace_all("[:space:]", "_")
colnames(clin) <- clin_colnames


# Tidy --------------------------------------------------------------------

clin_tidy <- clin %>% 
        dplyr::rename(id = "#patient_identifier",
               os_months = "overall_survival_(months)",
               dead = "overall_survival_status",
               response = "durable_clinical_benefit",
               recurred = "disease_free_status") %>%
        dplyr::select(-dosage, -`response_duration_(weeks)`, -cohort) %>% 
        mutate(dead = case_when(dead == "0:LIVING" ~ 0,
                                dead == "1:DECEASED" ~ 1,
                                T ~ NA_real_),
               recurred = case_when(recurred == "0:DiseaseFree" ~ 0,
                                recurred == "1:Recurred/Progressed" ~ 1,
                                T ~ NA_real_)) %>% 
        filter(id %in% colnames(rna)) %>% 
        column_to_rownames("id")


# Create Summarized Experiment --------------------------------------------

se <- SummarizedExperiment(rna, colData = DataFrame(clin_tidy))


# Get HGNC Symbol ---------------------------------------------------------

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ids <- getBM(attributes = c('hgnc_symbol','entrezgene_id', 'ensembl_gene_id'),
             filters = 'entrezgene_id',
             values = rownames(se),
             mart = mart)
ids2 <- ids[match(rownames(se), ids$entrezgene_id),]
mcols(se) <- cbind(mcols(se), ids2)


# Remove NA, blank HUGO symbols -------------------------------------------

se <- se[which(!is.na(rowData(se)$hgnc_symbol)),]
se <- se[which(rowData(se)$hgnc_symbol != ""),]
rowData(se)$hgnc_unique <- make.unique(rowData(se)$hgnc_symbol)
rownames(se) <- rowData(se)$hgnc_unique


# Logscale Counts ---------------------------------------------------------

assay(se, 2) <- apply(assay(se), 2, function(x) log2(x + 1))
assayNames(se) <- c("fpkm", "log_2")


# Write Data --------------------------------------------------------------

write_rds(se, "./data/van-allen/se.rds", compress = "gz")
