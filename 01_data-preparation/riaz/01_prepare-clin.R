
# Description -------------------------------------------------------------

# Prepare Riaz melanoma data
# Clinical data downloaded from here: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/
# And also in the 9th supplemental file of the paper

# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(readxl)


# Read in Data ------------------------------------------------------------

clin <- read_tsv("./data/riaz/GSE91061_series_matrix.txt", skip = 34, n_max = 45)
clin_2 <- read_excel("./data/riaz/NIHMS907788-supplement-9.xlsx", skip = 2)


# Tidy Annotations --------------------------------------------------------

tidy_annot <- clin %>%
        t() %>% 
        `colnames<-`(make.unique(.[1,])) %>% 
        as_tibble(rownames = "id") %>% 
        dplyr::slice(-1) %>% 
        mutate(patient = str_extract(id, "^[^_]*")) %>% 
        left_join(clin_2, by = c("patient" = "Patient"))

useful_info <- apply(tidy_annot, 2, function(x) sum(duplicated(x)))
useful_vals <- names(which(useful_info != 108))

tidy_annot <- tidy_annot %>% 
        dplyr::select(all_of(useful_vals)) %>% 
        dplyr::rename("geo" = `!Sample_geo_accession`, 
                      pre_or_on_treat = `!Sample_characteristics_ch1`,
                      response = `!Sample_characteristics_ch1.1`,
                      biosample = `!Sample_relation`,
                      sra = `!Sample_relation.1`) %>% 
        mutate(pre_or_on_treat = str_remove(pre_or_on_treat, "visit \\(pre or on treatment\\): "),
               pre_or_on_treat = as.factor(pre_or_on_treat),
               response = str_remove(response, "response: "),
               response = factor(response, levels= c("PRCR", "SD", "PD", "UNK")),
               sra = str_remove(sra, "SRA: "),
               biosample = str_remove(biosample, "BioSample: "))


# Write Tidy Clinical -----------------------------------------------------

write_rds(tidy_annot, "./data/riaz/tidy_clinical.Rds", compress = "gz")
