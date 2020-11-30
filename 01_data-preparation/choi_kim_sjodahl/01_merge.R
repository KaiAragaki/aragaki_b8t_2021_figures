
# Description -------------------------------------------------------------

# Merge Choi, Kim, and Sjodahl datasets


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(sva)
library(PCAtools)

# Read in Data ------------------------------------------------------------

kim_clin <- read_rds("./data/kim/01_tidy-clin.Rds")
kim_exp <- read_rds("./data/kim/01_tidy-expression.Rds")
kim_annot <- read_rds("./data/kim/01_tidy-annot.Rds")

sjodahl_clin <- read_rds("./data/sjodahl/01_tidy-clin.Rds")
sjodahl_exp <- read_rds("./data/sjodahl/01_tidy-expression.Rds")
sjodahl_annot <- read_rds("./data/sjodahl/01_tidy-annot.Rds")


choi_clin <- read_rds("./data/choi/01_tidy-clin.Rds")
choi_exp <- read_rds("./data/choi/01_tidy-expression.Rds")
choi_annot <- read_rds("./data/choi/01_tidy-annot.Rds")

annot <- full_join(kim_annot, sjodahl_annot, by = "nuID") %>% 
        full_join(choi_annot, by = "nuID")


# Select Needed Columns, Rows ---------------------------------------------

kim_clin <- kim_clin %>% 
        filter(`biological source:ch1` == "primary bladder cancer") %>% 
        filter(`invasiveness:ch1` == "muscle invasive") %>% 
        select(geo_accession, `AGE:ch1`, `cancer specific survival:ch1`, 
               `overall survival:ch1`, `SEX:ch1`, 
               `stage:ch1`, `survival month:ch1`)

kim_exp <- 
        kim_exp[ , which(colnames(kim_exp) %in% kim_clin$geo_accession)] %>% 
        as_tibble(rownames = "gene")

sjodahl_clin <- sjodahl_clin %>% 
        filter(!is.na(`time_to_dod_(months):ch1`)) %>% 
        filter(!(`tumor_stage:ch1` %in% c("Ta", "T1"))) %>% 
        select(geo_accession, `age:ch1`, `dod_event_(yes/no):ch1`, 
               `gender:ch1`, `tumor_stage:ch1`, `time_to_dod_(months):ch1`, `node_status:ch1`)

sjodahl_exp <- 
        sjodahl_exp[ , which(colnames(sjodahl_exp) %in% sjodahl_clin$geo_accession)] %>% 
        as_tibble(rownames = "gene")

choi_clin <- choi_clin %>% 
        filter(!is.na(`survival (mo):ch1`)) %>% 
        select(geo_accession, `age (at specimen collection):ch1`, 
               `dss censor:ch1`, `os censor:ch1`, Gender, `cstage:ch1`, 
               `survival (mo):ch1`)

choi_exp <- 
        choi_exp[ , which(colnames(choi_exp) %in% choi_clin$geo_accession)] %>% 
        as_tibble(rownames = "gene")


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
        
merged_exp <- full_join(kim_exp, sjodahl_exp, by = "gene") %>% 
        full_join(choi_exp, by = "gene")

merged_annot <- left_join(merged_exp, annot, by = c("gene" = "nuID")) %>% 
        mutate(symbol = case_when(is.na(symbol) & is.na(symbol.x) ~ NA_character_,
                                  is.na(symbol) & !is.na(symbol.x) ~ symbol.x,
                                  !is.na(symbol) & is.na(symbol.x) ~ symbol,
                                  symbol == symbol.x ~ symbol,
                                  symbol != symbol.x ~ paste(symbol, symbol.x, sep = "_"),
                                  T ~ NA_character_),
               symbol = case_when(is.na(symbol) & is.na(symbol.y) ~ NA_character_,
                                  is.na(symbol) & !is.na(symbol.y) ~ symbol.y,
                                  !is.na(symbol) & is.na(symbol.y) ~ symbol,
                                  symbol == symbol.y ~ symbol,
                                  symbol != symbol.y ~ paste(symbol, symbol.y, sep = "_"),
                                  T ~ NA_character_)) %>%
        select(-symbol.x, -symbol.y) %>% 
        relocate(contains("symbol"), .before = "gene")



signatures <- read_rds("./data/signatures/signatures.Rds")
signatures$b_cell[5] <- "KIAA0125" # Old HGNC used

query <- signatures %>% 
        unlist() %>% 
        unname() %>% 
        paste0(collapse = "$|")

sig_vec <- signatures %>% 
        unlist() %>% 
        unname()

probes_of_interest <- merged_annot[which(grepl(query, merged_annot$symbol)),] %>% 
        rowwise() %>% 
        mutate(num_na = sum(is.na(c_across(contains("GSM"))))) %>% 
        relocate(num_na, .before = symbol) %>% 
        group_by(symbol) %>% 
        arrange(symbol, num_na) %>% 
        fill(everything(), .direction = "downup") %>% 
        distinct(symbol, .keep_all = T)

# Merge probes with clin
tidy_probes <- t(probes_of_interest[-c(1:3)]) %>% 
        `colnames<-`(probes_of_interest$symbol) %>% 
        as_tibble(rownames = "geo_accession")

merged_long <- left_join(merged, tidy_probes, by = "geo_accession") %>% 
        pivot_longer(cols = BLK:TNFRSF17)

ggplot(merged_long, aes(x = value, fill = dataset)) + geom_histogram(alpha = 0.5) + facet_wrap(~name, scales = "free")


# Put the tidy probes back in to main expression matrix

probes <- probes_of_interest %>% 
        ungroup() %>% 
        select(-c(num_na, symbol))

expression <- filter(merged_annot, !(gene %in% probes_of_interest$gene)) %>% 
        bind_rows(probes_of_interest) %>% 
        group_by(symbol) %>% 
        arrange(symbol, num_na) %>% 
        fill(everything(), .direction = "downup") %>% 
        distinct(symbol, .keep_all = T)
        
        


expression_mat <- expression %>% 
        select(-gene, -num_na) %>% 
        column_to_rownames("symbol")

expression_mat_complete <- expression_mat[complete.cases(expression_mat),]

batch <- as.integer(as.factor(merged$dataset))

combat <- ComBat(as.matrix(expression_mat_complete), batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

combat_tib <- as_tibble(combat, rownames = "gene")

goi <- filter(combat_tib, gene %in% sig_vec) %>% 
        t() %>% 
        `colnames<-`(.[1,]) %>% 
        as_tibble(rownames = "sample") %>% 
        slice(-(1)) %>% 
        mutate(batch = batch) %>% 
        pivot_longer(BLK:TNFRSF17) %>% 
        mutate(value = as.numeric(value))

ggplot(goi, aes(value, fill = as.factor(batch))) + geom_density() + facet_wrap(~name, scales = "free") 

color_key <- data.frame(row.names = colnames(expression_mat_complete), group = merged$dataset)
biplot(pca(expression_mat_complete, metadata = color_key), colby = 'group', showLoadings = T, lab = NULL)
biplot(pca(combat, metadata = color_key), colby = 'group', showLoadings = T, lab = NULL)

# Write Merged Dataset -----------------------------------------------

write_tsv(merged, file = "./data/merged/kim_sjodahl_choi.tsv")
