
# Description -------------------------------------------------------------


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(DESeq2)
library(gt)
library(broom)


# Read in Data -----------------------------------------------------------------

blca <- read_rds("./data/tcga-blca/C_04_merge-gsva.rds")

clin <- as_tibble(colData(blca)) %>% 
        mutate(b.bin = if_else(b_cell > 0, "Hi", "Lo"),
               b.bin = factor(b.bin, levels = c("Lo", "Hi")),
               t.bin = if_else(cd8_rose > 0, "Hi", "Lo"),
               t.bin = factor(t.bin, levels = c("Lo", "Hi"))) %>% 
        unite(b8t, b.bin, t.bin, sep = ".", remove = F) %>% 
        mutate(sex = patient.gender,
               race = patient.race_list.race,
               b8t = factor(b8t, 
                            levels = c("Lo.Lo", "Lo.Hi", "Hi.Lo", "Hi.Hi"),
                            labels = c("Lo/Lo", "Lo/Hi", "Hi/Lo", "Hi/Hi"))) %>% 
        dplyr::select(sex, race, contains("metastat"), contains("stage"), 
               contains("age"), contains("smok"), contains("grade"), contains("lund"), 
               contains("mRNA"), "b8t", new_death, death_event, contains("sample.days_to_collection")) %>% 
        dplyr::select(!contains("file")) %>% 
        mutate(met_site = case_when(patient.metastatic_site_list.metastatic_site == "other" ~ patient.other_metastatic_site,
                                    T ~ patient.metastatic_site_list.metastatic_site)) %>% 
        unite(met_site, met_site, 
              patient.metastatic_site_list.metastatic_site.2, 
              patient.metastatic_site_list.metastatic_site.3, na.rm = T) %>% 
        mutate(met_site = if_else(met_site == "", NA_character_, met_site)) %>% 
        dplyr::select(-contains("prostate"), -contains("metastatic_site"), 
               -contains("metastatic_procedure"), -contains("gleason"), 
               -contains("dosage"), -contains("image"), -contains("omf")) %>% 
        dplyr::rename(stage = patient.stage_event.pathologic_stage,
                      stage_ct = patient.stage_event.tnm_categories.clinical_categories.clinical_t,
                      stage_m = patient.stage_event.tnm_categories.pathologic_categories.pathologic_m,
                      stage_n = patient.stage_event.tnm_categories.pathologic_categories.pathologic_n,
                      stage_pt = patient.stage_event.tnm_categories.pathologic_categories.pathologic_t,
                      age_smoking_start = patient.age_began_smoking_in_years,
                      grade = patient.neoplasm_histologic_grade,
                      age = patient.age_at_initial_pathologic_diagnosis,
                      days_to_collection = patient.samples.sample.days_to_collection,
                      pack_years = patient.number_pack_years_smoked) %>% 
        dplyr::select(-patient.stage_event.system_version, -stage) %>%  
        relocate(age, .after = age_smoking_start) %>% 
        mutate(stage_ct = case_when(str_detect(stage_ct, "t1") ~ "t1",
                                    str_detect(stage_ct, "t2") ~ "t2",
                                    str_detect(stage_ct, "t3") ~ "t3",
                                    str_detect(stage_ct, "t4") ~ "t4",
                                    T ~ stage_ct),
               stage_ct_bin = case_when(stage_ct == "t1" ~ "<T2",
                                        stage_ct %in% c("t2", "t3", "t4") ~ ">=T2",
                                        T ~ stage_ct),
               stage_pt = case_when(str_detect(stage_pt, "t1") ~ "t1",
                                    str_detect(stage_pt, "t2") ~ "t2",
                                    str_detect(stage_pt, "t3") ~ "t3",
                                    str_detect(stage_pt, "t4") ~ "t4",
                                    T ~ stage_pt),
               met_site_bin = case_when(met_site == "lymph node only" ~ "Lymph Node",
                                        met_site == "none" ~ "None",
                                        is.na(met_site) ~ NA_character_,
                                        T ~ "Other"),
               met_site_bin = factor(met_site_bin, levels = c("None", "Lymph Node", "Other")),
               node_bin = case_when(stage_n %in% c("n1", "n2", "n3") ~ ">N0",
                                    stage_n == "n0" ~ "N0",
                                    stage_n == "nx" ~ "Nx",
                                    T ~ stage_n),
               node_bin = factor(node_bin, levels = c("N0", ">N0", "Nx")),
               stage_m = case_when(stage_m == "m0" ~ "M0",
                                   stage_m == "m1" ~ "M1",
                                   stage_m == "mx" ~ "Mx",
                                   T ~ stage_m),
               mRNA.cluster = factor(mRNA.cluster, 
                                     levels = c("Basal_squamous", "Luminal", "Luminal_infiltrated", "Luminal_papillary", "Neuronal"),
                                     labels = c("Basal Squamous", "Luminal", "Luminal Infiltrated", "Luminal Papillary", "Neuronal")),
               grade = factor(grade, levels = c("low grade", "high grade"), labels = c("Low", "High")),
               sex = factor(sex,
                            levels = c("male", "female"),
                            labels = c("Male", "Female"))) %>% 
        select(-stage_ct, -stage_pt, -met_site, -stage_n) %>% 
        mutate(across(where(is.character), as.factor))

# Didn't do stage_pt_bin because it was causing convergence issues due to low <t2.

# Univariable ------------------------------------------------------------------

## Log-rank Test ---------------------------------------------------------------

# Cutoff for consideration in multivariable analysis (denoted 'sig'): 
# p < 0.15

lr <- function(var) {
        name <- var
        var <- clin[[var]]
        coxph(Surv(new_death, death_event) ~ var, data = clin) %>% 
                glance() %>% 
                mutate(name = name) %>% 
                relocate(name)
}

vars <- clin %>% 
        dplyr::select(-new_death, -death_event) %>% 
        colnames()

df <- data.frame()
for (i in seq_along(vars)) {
        df <- rbind(df, as.vector(lr(vars[i])))
}

df <- df %>% 
        select(name, p.value.sc) %>% 
        mutate(stars = case_when(p.value.sc < 0.001 ~ "(***)",
                                 p.value.sc < 0.01 ~ "(**)",
                                 p.value.sc < 0.05 ~ "(*)",
                                 p.value.sc < 0.15 ~ "(#)",
                                 T ~ "(NS)"))

df_sig <- df %>% 
        filter(p.value.sc < 0.15 | name %in% c("sex", "b8t")) %>% 
        unite(term_2, name, stars, remove = F, sep = " ") %>% 
        select(-stars)


# Cox with Levels ---------------------------------------------------------

# Now look at coxph ratios for each level that reached ss

wl <- function(var) {
        feature <- var
        var <- clin[[var]]
        coxph(Surv(new_death, death_event) ~ var, data = clin) %>% 
                tidy(exponentiate = TRUE, conf.int = TRUE) %>%  
                mutate(feature = feature) %>% 
                relocate(feature) %>% 
                retidy()
}

colnames_search <- paste(colnames(clin), collapse = "|")

make_reflevel_table <- function(df) {
        feature <- unique({{df$feature}})
        if(is.factor(clin[[feature]])) {
                ref_df <- data.frame(feature = feature, 
                                     term = levels(clin[[feature]])[1], 
                                     estimate = 1, std.error = NA, p.value = NA)
                df <- bind_rows(ref_df, df)
        }
        df
}

retidy <- function(df) {
        df %>% 
                mutate(term = str_remove(term, "^var")) %>% 
                select(-statistic) %>% 
                make_reflevel_table()
}

vars <- df_sig$name
df <- data.frame()
for (i in seq_along(vars)) {
        df <- rbind(df, as.vector(wl(vars[i])))
}


# Make Table ---------------------------------------------------------

df %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "***",
                                 p.value < 0.01 ~ "**",
                                 p.value < 0.05 ~ "*",
                                 is.na(p.value) ~ NA_character_,
                                 T ~ "NS")) %>% 
        left_join(df_sig, by = c("feature" = "name")) %>%
        dplyr::rename("HR" = estimate,
                      "p-value" = p.value,
                      " " = stars) %>%
        select(-feature, -p.value.sc, -std.error) %>% 
        relocate("p-value", .before =  " ") %>% 
        mutate(term_2 = str_replace(term_2, "stage_m", "Metastasis"),
               term_2 = str_replace(term_2, "age", "Age"),
               term_2 = str_replace(term_2, "grade", "Grade"),
               term_2 = str_replace(term_2, "mRNA.cluster", "TCGA Subtype"),
               term_2 = str_replace(term_2, "days_to_collection", "Days to Collection"),
               term_2 = str_replace(term_2, "met_site_bin", "Site of Metastsis"),
               term_2 = str_replace(term_2, "node_bin", "Node Status"),
               term_2 = str_replace(term_2, "sex", "Sex"),
               term_2 = str_replace(term_2, "b8t", "B8T")) %>%
        gt(groupname_col = "term_2",
           rowname_col = "term") %>% 
        fmt_number(columns = 2:4, drop_trailing_zeros = T) %>% 
        fmt_missing(columns = 1:6, missing_text = "") %>%
        fmt_scientific(columns = 5, decimals = 2) %>% 
        cols_merge(columns = c("conf.low", "conf.high"),
                   pattern = "{1}-{2}") %>%
        cols_label(conf.low = "CI (95%)") %>%
        cols_align("center") %>% 
        tab_style(style = cell_text(align = "right"), 
                  locations = cells_stub()) %>% 
        tab_header(title = "Univariable Hazard Ratios") %>% 
        tab_footnote("Features not trending significant (P < 0.15) by logrank test: Sex, Race, Age Started Smoking, B8T, Clinical Stage",
                     locations = cells_column_labels("p-value")) %>% 
        tab_footnote("p-value by logrank test (*** < 0.001, ** < 0.01, * < 0.05, # < 0.15)",
                     location = cells_row_groups(starts_with("stage_m"))) %>% 
        cols_width(vars(`p-value`, `conf.low`) ~ px(130),
                   vars(HR, " ") ~ px(60)) %>%
        tab_options(table.width = px(600)) %>% 
        gtsave("./figures/tables/risk_table_uni_blca.png") 


# Multivariate -------------------------------------------------------

complete_clin <- clin %>% 
        select(b8t, sex, stage_m, age, grade, mRNA.cluster, days_to_collection, met_site_bin, node_bin, new_death, death_event)

complete_clin <- complete_clin[complete.cases(complete_clin),]

retidy <- function(df) {
        df %>% 
                mutate(feature = str_extract(term, colnames_search),
                       term = str_remove(term, colnames_search)) %>% 
                group_by(feature) %>% 
                select(-statistic)
}

make_reflevel_table <- function(df) {
        feature <- unique({{df$feature}})
        if(is.factor(complete_clin[[feature]])) {
                ref_df <- data.frame(feature = feature, 
                                     term = levels(complete_clin[[feature]])[1], 
                                     estimate = 1, std.error = NA, p.value = NA)
                df <- bind_rows(ref_df, df)
        }
        df
}

mv <- coxph(Surv(new_death, death_event) ~ b8t*sex + stage_m + age + grade + mRNA.cluster + days_to_collection + met_site_bin + node_bin, data = complete_clin) %>% 
        tidy(exponentiate = T) %>%
        retidy() %>% 
        group_by(feature) %>% 
        group_split() %>% 
        lapply(make_reflevel_table) %>% 
        bind_rows() %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "(***)",
                                 p.value < 0.01 ~ "(**)",
                                 p.value < 0.05 ~ "(*)",
                                 T ~ "(NS)")) %>% 
        filter(!(feature %in% c("grade", "met_site_bin", "node_bin", "stage_m")))

# Not significant:
# Grade, Met_Site_Bin, Node_Bin, Sex, Stage_M


# Make Table ---------------------------------------------------------

mv %>% 
        dplyr::rename("HR" = estimate,
                      "SE" = std.error,
                      "p-value" = p.value,
                      " " = stars) %>%
        mutate(feature = str_replace(feature, "b8t", "B8T"),
               term = str_replace(term, "Hi.Lo", "Hi/Lo"),
               term = str_replace(term, "Lo.Lo", "Lo/Lo"),
               term = str_replace(term, "Lo.Hi", "Lo/Hi"),
               term = str_replace(term, "Hi.Hi", "Hi/Hi"),
               term = str_replace(term, "(?<!fe)male", "Male"),
               term = str_replace(term, "female", "Female"),
               term = str_remove(term, "sex"),
               feature = case_when(str_detect(feature, "age") ~ "Age",
                                   str_detect(feature, "mRNA.cluster") ~ "TCGA Subtype",
                                   str_detect(feature, "days_to_collection") ~ "Days to Collection",
                                   str_detect(feature, "sex") ~ "Sex",
                                   T ~ feature)) %>% 
        gt(groupname_col = "feature",
           rowname_col = "term") %>% 
        fmt_number(columns = 2:4, drop_trailing_zeros = T) %>% 
        fmt_scientific(columns = 5, decimals = 2) %>% 
        cols_align("center") %>% 
        tab_style(style = cell_text(align = "right"), 
                  locations = cells_stub()) %>% 
        fmt_missing(columns = 1:6, missing_text = "") %>%
        tab_header(title = "Multivariate Hazard Ratios") %>% 
        tab_footnote("B8T and Sex are not significant (p < 0.15) in the univariable analysis, but are included here to be complete",
                     locations = cells_column_labels("p-value")) %>% 
        cols_width(vars(`p-value`) ~ px(130),
                   vars(HR, SE, " ") ~ px(60)) %>%
        tab_options(table.width = px(500)) %>% 
        gtsave("./figures/tables/risk_table_multi_blca.png") 

