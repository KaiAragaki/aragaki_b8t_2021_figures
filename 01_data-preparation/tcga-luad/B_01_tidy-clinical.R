
# Description -------------------------------------------------------------

# This file tidies the TCGA LUAD clinical data


# Prepare Workspace -------------------------------------------------------

library(tidyverse)


# Read in Data ------------------------------------------------------------

clinical <- download.file("http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.Merge_Clinical.Level_1.2016012800.0.0.tar.gz", destfile = "./data/tcga-luad/luad-clin.tar.gz")
untar("./data/tcga-luad/luad-clin.tar.gz", exdir = "./data/tcga-luad/")

clinical <- read_tsv("./data/tcga-luad/gdac.broadinstitute.org_LUAD.Merge_Clinical.Level_1.2016012800.0.0/LUAD.clin.merged.txt", col_names = F) %>% 
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
                dplyr::select(collapsed)
}

death <- 
        clin_filtered %>%
        dplyr::select(contains("days_to_death")) %>% 
        get_latest()

fl <- 
        clin_filtered %>%
        dplyr::select(contains("days_to_last_followup")) %>% 
        get_latest()

all_clin <- tibble(death, fl, .name_repair = "unique") %>% 
        `colnames<-`(c('death_days', 'followUp_days')) %>% 
        mutate(new_death = if_else(is.na(death_days), followUp_days, death_days),
               death_event = if_else(clinical$patient.vital_status == 'dead'|!is.na(death_days), 1, 0)) %>% 
        cbind(clin_filtered)


# Save File ---------------------------------------------------------------

write_csv(all_clin, "./data/tcga-luad/B_01_tidy-clinical.csv")

