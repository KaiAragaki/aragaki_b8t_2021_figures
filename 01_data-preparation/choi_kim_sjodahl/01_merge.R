
# Description -------------------------------------------------------------

# Merge Choi, Kim, and Sjodahl datasets


# Prepare Workspace -------------------------------------------------------

library(tidyverse)


# Read in Data ------------------------------------------------------------

kim_clin <- read_rds("./data/kim/01_tidy-clin.Rds")
kim_exp <- read_rds("./data/kim/01_tidy-expression.Rds")

sjodahl_clin <- read_rds("./data/sjodahl/01_tidy-clin.Rds")
sjodahl_exp <- read_rds("./data/sjodahl/01_tidy-expression.Rds")

choi_clin <- read_rds("./data/choi/01_tidy-clin.Rds")
choi_exp <- read_rds("./data/choi/01_tidy-expression.Rds")


# Select Needed Columns, Rows ---------------------------------------------

kim_clin <- kim_clin %>% 
        filter(`biological source:ch1` == "primary bladder cancer") %>% 
        filter(`invasiveness:ch1` == "muscle invasive") %>% 
        select(geo_accession, `AGE:ch1`, `cancer specific survival:ch1`, 
               `overall survival:ch1`, `SEX:ch1`, 
               `stage:ch1`, `survival month:ch1`)

kim_exp <- 
        kim_exp[ , which(colnames(kim_exp) %in% kim_clin$geo_accession)]

sjodahl_clin <- sjodahl_clin %>% 
        filter(!is.na(`time_to_dod_(months):ch1`)) %>% 
        filter(!(`tumor_stage:ch1` %in% c("Ta", "T1"))) %>% 
        select(geo_accession, `age:ch1`, `dod_event_(yes/no):ch1`, 
               `gender:ch1`, `tumor_stage:ch1`, `time_to_dod_(months):ch1`, `node_status:ch1`)

sjodahl_exp <- 
        sjodahl_exp[ , which(colnames(sjodahl_exp) %in% sjodahl_clin$geo_accession)]

choi_clin <- choi_clin %>% 
        filter(!is.na(`survival (mo):ch1`)) %>% 
        select(geo_accession, `age (at specimen collection):ch1`, 
               `dss censor:ch1`, `os censor:ch1`, Gender, `cstage:ch1`, 
               `survival (mo):ch1`)

choi_exp <- 
        choi_exp[ , which(colnames(choi_exp) %in% choi_clin$geo_accession)]


# Rename Columns to Match -------------------------------------------------

kim_clin <- rename(kim_clin,
                   age = "AGE:ch1",
                   censor_dss = "cancer specific survival:ch1",
                   censor_os = "overall survival:ch1",
                   sex = "SEX:ch1",
                   stage = "stage:ch1",
                   survival_month = "survival month:ch1")

sjodahl_clin <- sjodahl_clin %>% 
        mutate(censor_dss = NA) %>% 
        rename(age = "age:ch1",
               censor_os = "dod_event_(yes/no):ch1",
               sex = "gender:ch1",
               survival_month = "time_to_dod_(months):ch1") %>%
        unite(stage, `tumor_stage:ch1`, `node_status:ch1`, sep = "", na.rm = T) %>% 
        relocate(censor_dss, .before = censor_os)

choi_clin <- rename(choi_clin,
                    age = "age (at specimen collection):ch1",
                    censor_dss = "dss censor:ch1",
                    censor_os = "os censor:ch1",
                    sex = "Gender",
                    stage = "cstage:ch1",
                    survival_month = "survival (mo):ch1")


# Homogenize Factor Levels and Labels --------------------------------

kim_clin <- kim_clin %>% 
        mutate(censor_dss = if_else(censor_dss == "survival", 0, 1), # Checked to make sure it's just "survival" and "death" with no NAs
               censor_os = if_else(censor_os == "survival", 0, 1),
               dataset = "kim") 

sjodahl_clin <- sjodahl_clin %>% 
        mutate(censor_os = if_else(censor_os == "no", 0, 1),
               dataset = "sjodahl")

choi_clin <- choi_clin %>% 
        mutate(censor_dss = if_else(censor_dss == "uncensored", 0, 1),
               censor_os = if_else(censor_os == "uncensored", 0, 1),
               sex = if_else(tolower(sex) == "male", "M", "F"),
               dataset = "choi") 


# Merge Data ---------------------------------------------------------

merged <- bind_rows(kim_clin, sjodahl_clin, choi_clin) %>% 
        mutate(stage = str_replace(stage, "^cT", "T"),
               stage = str_replace(stage, "pN", "N"),
               stage = str_replace(stage, "^pT", "T"),
               stage = str_remove(stage, "M.*$"),
               stage = str_replace(stage, "N", "_N")) %>%
        separate(stage, into = c("tumor", "node"), sep = "_")
        

# Write Merged Dataset -----------------------------------------------

write_tsv(merged, file = "./data/merged/kim_sjodahl_choi.tsv")
