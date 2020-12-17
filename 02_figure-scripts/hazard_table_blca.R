
# Description -------------------------------------------------------------


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(DESeq2)
library(gt)
library(broom)


# Read in Data ------------------------------------------------------------

blca <- read_rds("./data/tcga-blca/C_04_merge-gsva.rds")

clin <- as_tibble(colData(blca)) %>% 
        mutate(b.bin = if_else(b_cell > 0, "hi", "lo"),
               t.bin = if_else(cd8_rose > 0, "hi", "lo")) %>% 
        unite(b8t, b.bin, t.bin, sep = ".", remove = F) %>% 
        mutate(sex = patient.gender,
               race = patient.race_list.race) %>% 
        select(sex, race, contains("metastat"), contains("stage"), 
               contains("age"), contains("grade"), contains("lund"), 
               contains("mRNA"), "b8t", new_death, death_event, contains("sample.days_to_collection")) %>% 
        select(!contains("file")) %>% 
        mutate(met_site = case_when(patient.metastatic_site_list.metastatic_site == "other" ~ patient.other_metastatic_site,
                                    T ~ patient.metastatic_site_list.metastatic_site)) %>% 
        unite(met_site, met_site, 
              patient.metastatic_site_list.metastatic_site.2, 
              patient.metastatic_site_list.metastatic_site.3, na.rm = T) %>% 
        select(-contains("prostate"), -contains("metastatic_site"), 
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
                      days_to_collection = patient.samples.sample.days_to_collection) %>% 
        select(-patient.stage_event.system_version) %>% 
        relocate(stage, .after = stage_pt) %>% 
        relocate(age, .after = age_smoking_start) %>% 
        mutate(stage_ct = case_when(str_detect(stage_ct, "t1") ~ "t1",
                                    str_detect(stage_ct, "t2") ~ "t2",
                                    str_detect(stage_ct, "t3") ~ "t3",
                                    str_detect(stage_ct, "t4") ~ "t4",
                                    T ~ stage_ct),
               stage_pt = case_when(str_detect(stage_pt, "t1") ~ "t1",
                                    str_detect(stage_pt, "t2") ~ "t2",
                                    str_detect(stage_pt, "t3") ~ "t3",
                                    str_detect(stage_pt, "t4") ~ "t4",
                                    T ~ stage_pt)) %>% 
        mutate(across(where(is.character), as.factor))


# Univariate ---------------------------------------------------------


# Anovas -------------------------------------------------------------



# Cutoff for consideration in multivariate analysis ("denoted 'sig'"): 
# p < 0.15

a <- coxph(Surv(new_death, death_event) ~ sex, data = clin) %>% 
        anova() %>%  
        tidy()
b <- coxph(Surv(new_death, death_event) ~ race, data = clin) %>% 
        anova() %>% 
        tidy()
c <- coxph(Surv(new_death, death_event) ~ age_smoking_start, data = clin) %>% 
        anova() %>% 
        tidy()
d <- coxph(Surv(new_death, death_event) ~ mRNA.cluster, data = clin) %>% 
        anova() %>% 
        tidy()
# e <- coxph(Surv(new_death, death_event) ~ stage, data = clin) %>% 
#         anova() %>% 
#         tidy()

clin_2 <- filter(clin, stage != "stage i")

e.2 <- coxph(Surv(new_death, death_event) ~ stage, data = clin_2) %>% 
        anova() %>% 
        tidy()

f <- coxph(Surv(new_death, death_event) ~ stage_ct, data = clin) %>% 
        anova() %>% 
        tidy()

g <- coxph(Surv(new_death, death_event) ~ stage_m, data = clin) %>% 
        anova() %>% 
        tidy()
h <- coxph(Surv(new_death, death_event) ~ stage_n, data = clin) %>% 
        anova() %>% 
        tidy()
# i <- coxph(Surv(new_death, death_event) ~ stage_pt, data = clin) %>% 
#         anova() %>% 
#         tidy()

pt_table <- table(clin$stage_pt)
few <- names(pt_table[which(pt_table > 10)])

clin_3 <- filter(clin, stage_pt %in% few)

i.2 <- coxph(Surv(new_death, death_event) ~ stage_pt, data = clin_3) %>% 
        anova() %>% 
        tidy()

j <- coxph(Surv(new_death, death_event) ~ age, data = clin) %>% 
        anova() %>% 
        tidy()
k <- coxph(Surv(new_death, death_event) ~ grade, data = clin) %>% 
        anova() %>% 
        tidy()
l <- coxph(Surv(new_death, death_event) ~ b8t, data = clin) %>% 
        anova() %>% 
        tidy()
m <- coxph(Surv(new_death, death_event) ~ days_to_collection, data = clin) %>% 
        anova() %>% 
        tidy()

signif_uni <- 
        bind_rows(a, b, c, d, e.2, f, g, h, i.2, j, k, l, ) %>% 
        filter(term != "NULL") %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "(***)",
                                 p.value < 0.01 ~ "(**)",
                                 p.value < 0.05 ~ "(*)",
                                 p.value < 0.15 ~ "(#)",
                                 T ~ "(NS)")) %>% 
        filter(p.value < 0.15) %>% 
        unite(term_2, term, stars, sep = " ", remove = F)
        

# Cox with Levels ----------------------------------------------------

# Now look at coxph ratios for each level that reached ss

colnames_search <- paste(colnames(clin), collapse = "|")

make_reflevel_table <- function(df) {
        name <- unique({{df$feature}})
        if(is.factor(clin[[name]])) {
                ref_df <- data.frame(feature = name, term = levels(clin[[name]])[1], estimate = 1, std.error = NA, p.value = NA)
                df <- bind_rows(ref_df, df)
        }
        df
}

retidy <- function(df) {
        df %>% 
                mutate(feature  = str_extract(term, colnames_search),
                       term = str_remove(term, colnames_search)) %>% 
                relocate(feature, .before = everything()) %>% 
                select(-statistic) %>% 
                make_reflevel_table()
}

a.1 <- coxph(Surv(new_death, death_event) ~ mRNA.cluster, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
b.1 <- coxph(Surv(new_death, death_event)~ stage, data = clin_2) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
c.1 <- coxph(Surv(new_death, death_event)~ stage_ct, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
d.1 <- coxph(Surv(new_death, death_event)~ stage_m, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
e.1 <- coxph(Surv(new_death, death_event)~ stage_n, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
f.1 <- coxph(Surv(new_death, death_event)~ stage_pt, data = clin_3) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
g.1 <- coxph(Surv(new_death, death_event)~ age, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
h.1 <- coxph(Surv(new_death, death_event)~ grade, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
i.1 <- coxph(Surv(new_death, death_event)~ b8t, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
j.1 <- coxph(Surv(new_death, death_event)~ days_to_collection, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()

# Did not reach SS (p < 0.15):
# gender
# race
# age_smoking_start


# Make Table ---------------------------------------------------------

bind_rows(a.1, b.1, c.1, d.1, e.1, f.1, g.1, h.1, i.1, j.1) %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "***",
                                 p.value < 0.01 ~ "**",
                                 p.value < 0.05 ~ "*",
                                 is.na(p.value) ~ NA_character_,
                                 T ~ "NS"))%>% 
        dplyr::rename("HR" = estimate,
                      "SE" = std.error,
                      "p-value" = p.value,
                      " " = stars) %>%
        left_join(signif_uni, by = c("feature" = "term")) %>%
        select(-(logLik:stars), -feature) %>% 
        mutate(term_2 = str_replace_all(term_2, "\\.|_", " ")) %>% 
        gt(groupname_col = "term_2",
           rowname_col = "term") %>% 
        fmt_number(columns = 2:3) %>% 
        fmt_scientific(columns = 4, decimals = 2) %>% 
        cols_align("center") %>% 
        tab_style(style = cell_text(align = "right"), 
                  locations = cells_stub()) %>% 
        fmt_missing(columns = 1:6, missing_text = "") %>%
        tab_header(title = "Univariate Hazard Ratios") %>% 
        tab_footnote("Features not trending significant (P < 0.15) by ANOVA: Sex, Race, Age-Smoking-Start",
                     locations = cells_column_labels("p-value")) %>% 
        tab_footnote("p-value by ANOVA (*** < 0.001, ** < 0.01, * < 0.05, # < 0.15)",
                     location = cells_row_groups(starts_with("mRNA"))) %>% 
        tab_footnote("Missing values are due to low sample number/unreliable estimates",
                     location = cells_row_groups(matches("^stage \\(|^stage pt"))) %>% 
        cols_width(vars(`p-value`) ~ px(130),
                   vars(HR, SE, " ") ~ px(60)) %>%
        tab_options(table.width = px(500)) %>% 
        gtsave("./figures/tables/blca_risk_table_uni.png") 


# Multivariate -------------------------------------------------------

# Adjusting for age

# Anovas -------------------------------------------------------------


a <- coxph(Surv(new_death, death_event) ~ age + sex + b8t, data = clin) %>% 
        anova() %>%  
        tidy()

b <- coxph(Surv(new_death, death_event) ~ age + mRNA.cluster + b8t, data = clin) %>% 
        anova() %>% 
        tidy()

c <- coxph(Surv(new_death, death_event) ~ age + stage + b8t, data = clin_2) %>% 
        anova() %>% 
        tidy()

d <- coxph(Surv(new_death, death_event) ~ age + stage_ct + b8t, data = clin) %>% 
        anova() %>% 
        tidy()

e <- coxph(Surv(new_death, death_event) ~ age + stage_m + b8t, data = clin) %>% 
        anova() %>% 
        tidy()
f <- coxph(Surv(new_death, death_event) ~ age + stage_n + b8t, data = clin) %>% 
        anova() %>% 
        tidy()

g <- coxph(Surv(new_death, death_event) ~ age + stage_pt + b8t, data = clin_3) %>% 
        anova() %>% 
        tidy()

h <- coxph(Surv(new_death, death_event) ~ age + grade + b8t, data = clin) %>% 
        anova() %>% 
        tidy()

i <- coxph(Surv(new_death, death_event) ~ age + days_to_collection + b8t, data = clin) %>% 
        anova() %>% 
        tidy()

signif_multi <- 
        bind_rows(a, b, c, d, e, f, g, h, i) %>% 
        filter(term != "NULL") %>% 
        filter(term != "age") %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "(***)",
                                 p.value < 0.01 ~ "(**)",
                                 p.value < 0.05 ~ "(*)",
                                 T ~ "(NS)")) %>% 
        unite(term_2, term, stars, sep = " ", remove = F)

# In all instances, b8t fails to reach significance when placed as the first predictor
# However, when adjusting for mRNA.cluster, b8t shows significance




# Cox with Levels ----------------------------------------------------

# Now look at coxph ratios for each level that reached ss

colnames_search <- paste(colnames(clin), collapse = "|")

make_reflevel_table <- function(df) {
        name <- unique({{df$feature}})
        if(is.factor(clin[[name]])) {
                ref_df <- data.frame(feature = name, term = levels(clin[[name]])[1], estimate = 1, std.error = NA, p.value = NA)
                df <- bind_rows(ref_df, df)
        }
        df
}

retidy <- function(df) {
        df %>% 
                mutate(feature  = str_extract(term, colnames_search),
                       term = str_remove(term, colnames_search)) %>% 
                relocate(feature, .before = everything()) %>% 
                select(-statistic) %>% 
                make_reflevel_table()
}

a.1 <- coxph(Surv(new_death, death_event) ~ mRNA.cluster, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
b.1 <- coxph(Surv(new_death, death_event)~ stage, data = clin_2) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
c.1 <- coxph(Surv(new_death, death_event)~ stage_ct, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
d.1 <- coxph(Surv(new_death, death_event)~ stage_m, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
e.1 <- coxph(Surv(new_death, death_event)~ stage_n, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
f.1 <- coxph(Surv(new_death, death_event)~ stage_pt, data = clin_3) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
g.1 <- coxph(Surv(new_death, death_event)~ age, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
h.1 <- coxph(Surv(new_death, death_event)~ grade, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
i.1 <- coxph(Surv(new_death, death_event)~ b8t, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
j.1 <- coxph(Surv(new_death, death_event)~ days_to_collection, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()

# Did not reach SS (p < 0.15):
# gender
# race
# age_smoking_start


# Make Table ---------------------------------------------------------

bind_rows(a.1, b.1, c.1, d.1, e.1, f.1, g.1, h.1, i.1, j.1) %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "***",
                                 p.value < 0.01 ~ "**",
                                 p.value < 0.05 ~ "*",
                                 is.na(p.value) ~ NA_character_,
                                 T ~ "NS"))%>% 
        dplyr::rename("HR" = estimate,
                      "SE" = std.error,
                      "p-value" = p.value,
                      " " = stars) %>%
        left_join(signif_uni, by = c("feature" = "term")) %>%
        select(-(logLik:stars), -feature) %>% 
        mutate(term_2 = str_replace_all(term_2, "\\.|_", " ")) %>% 
        gt(groupname_col = "term_2",
           rowname_col = "term") %>% 
        fmt_number(columns = 2:3) %>% 
        fmt_scientific(columns = 4, decimals = 2) %>% 
        cols_align("center") %>% 
        tab_style(style = cell_text(align = "right"), 
                  locations = cells_stub()) %>% 
        fmt_missing(columns = 1:6, missing_text = "") %>%
        tab_header(title = "Univariate Hazard Ratios") %>% 
        tab_footnote("Features not trending significant (P < 0.15) by ANOVA: Sex, Race, Age-Smoking-Start",
                     locations = cells_column_labels("p-value")) %>% 
        tab_footnote("p-value by ANOVA (*** < 0.001, ** < 0.01, * < 0.05, # < 0.15)",
                     location = cells_row_groups(starts_with("mRNA"))) %>% 
        tab_footnote("Missing values are due to low sample number/unreliable estimates",
                     location = cells_row_groups(matches("^stage \\(|^stage pt"))) %>% 
        cols_width(vars(`p-value`) ~ px(130),
                   vars(HR, SE, " ") ~ px(60)) %>%
        tab_options(table.width = px(500)) %>% 
        gtsave("./figures/tables/blca_risk_table_uni.png") 






# Forward Selection --------------------------------------------------


# Fwd selection, B8T as predictor, not forcing first -----------------

a <- coxph(Surv(new_death, death_event)~mRNA.cluster, clin)
b <- coxph(Surv(new_death, death_event)~stage, clin_2)
c <- coxph(Surv(new_death, death_event)~stage_ct, clin)
d <- coxph(Surv(new_death, death_event)~stage_m, clin)
e <- coxph(Surv(new_death, death_event)~stage_n, clin)
f <- coxph(Surv(new_death, death_event)~stage_pt, clin_3)
g <- coxph(Surv(new_death, death_event)~age, clin)
h <- coxph(Surv(new_death, death_event)~grade, clin)
i <- coxph(Surv(new_death, death_event)~b8t, clin)
j <- coxph(Surv(new_death, death_event)~days_to_collection, clin)

anova(a) # 6.98e-4
anova(b) # 1.50e-7
anova(c) # 0.13
anova(d) # 4.81e-3
anova(e) # 9.38e-6
anova(f) # 6.93e-5
anova(g) # 2.00e-5
anova(h) # 0.07
anova(i) # 0.14
anova(j) # 1.90e-2


a.1 <- coxph(Surv(new_death, death_event) ~ stage + mRNA.cluster, clin_2)
b.1 <- coxph(Surv(new_death, death_event) ~ stage + stage_ct, clin_2)
c.1 <- coxph(Surv(new_death, death_event) ~ stage + stage_m, clin_2)
d.1 <- coxph(Surv(new_death, death_event) ~ stage + stage_n, clin_2)
e.1 <- coxph(Surv(new_death, death_event) ~ stage + stage_pt, clin_3)
f.1 <- coxph(Surv(new_death, death_event) ~ stage + age, clin_2)
g.1 <- coxph(Surv(new_death, death_event) ~ stage + grade, clin_2)
h.1 <- coxph(Surv(new_death, death_event) ~ stage + b8t, clin_2)
i.1 <- coxph(Surv(new_death, death_event) ~ stage + days_to_collection, clin_2)

anova(a.1) # 0.054
anova(b.1) # 0.156
anova(c.1) # 0.122
anova(d.1) # 0.775
anova(e.1) # 0.243
anova(f.1) # 1.25e-4
anova(g.1) # 0.378
anova(h.1) # 0.381
anova(i.1) # 4.47e-3
