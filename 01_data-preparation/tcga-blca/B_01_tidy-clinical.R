
# Description -------------------------------------------------------------

# This file tidies the TCGA BLCA clinical data
# Last downloaded from Broad Firehose 2020-10-06
# http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BLCA/20160128/gdac.broadinstitute.org_BLCA.Merge_Clinical.Level_1.2016012800.0.0.tar.gz

# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)


# Read in Data ------------------------------------------------------------

download.file("http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BLCA/20160128/gdac.broadinstitute.org_BLCA.Merge_Clinical.Level_1.2016012800.0.0.tar.gz", 
              destfile = "./data/tcga-blca/gdac.broadinstitute.org_BLCA.Merge_Clinical.Level_1.2016012800.0.0.tar.gz")

untar("./data/tcga-blca/gdac.broadinstitute.org_BLCA.Merge_Clinical.Level_1.2016012800.0.0.tar.gz", exdir = "./data/tcga-blca/")

clinical <- read_tsv("./data/tcga-blca/gdac.broadinstitute.org_BLCA.Merge_Clinical.Level_1.2016012800.0.0/BLCA.clin.merged.txt", col_names = F) %>% 
        t()

colnames(clinical) <- clinical[1,]

clinical <- clinical %>% 
        as_tibble()

clinical <- clinical[-1,]


# Remove Columns with Only One Unique Entry -------------------------------

is_just_one <- function(x) {
        length(unique(x)) == 1
}

res <- apply(clinical, 2, is_just_one)

clin_filtered <- clinical[,!res] %>% 
        as_tibble()


# Create Survival Time Column ---------------------------------------------

get_latest <- function(data) {
        rowwise(data) %>% 
                mutate(across(everything(), as.numeric)) %>% 
                mutate(collapsed = if_else(!is.infinite(max(c_across(), na.rm = T)), 
                                           max(c_across(), na.rm = T), 
                                           NA_real_)) %>% 
                select(collapsed)
}

death <- 
        clin_filtered %>%
        select(contains("days_to_death")) %>% 
        get_latest()

fl <- 
        clin_filtered %>%
        select(contains("days_to_last_followup")) %>% 
        get_latest()

all_clin <- tibble(death, fl, .name_repair = "unique") %>% 
        `colnames<-`(c('death_days', 'followUp_days')) %>% 
        mutate(new_death = if_else(is.na(death_days), followUp_days, death_days),
               death_event = if_else(clinical$patient.vital_status == 'dead'|!is.na(death_days), 1, 0)) %>% 
        cbind(clin_filtered)


# Save File ---------------------------------------------------------------

write_csv(all_clin, "./data/tcga-blca/B_01_tidy-clinical.csv")
