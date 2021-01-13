
# Description -------------------------------------------------------------


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(DESeq2)
library(gt)
library(broom)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")

clin <- as_tibble(colData(imvigor)) %>% 
        mutate(b.bin = if_else(b_cell > 0, "hi", "lo"),
               b.bin = factor(b.bin, levels = c("lo", "hi")),
               t.bin = if_else(cd8_rose > 0, "hi", "lo"),
               t.bin = factor(t.bin, levels = c("lo", "hi")),) %>% 
        unite(b8t, b.bin, t.bin, sep = ".", remove = F) %>% 
        mutate(gender = factor(gender, levels = c("male", "female")),
               IC.Level = factor(IC.Level, levels = c("IC0", "IC1", "IC2+")),
               Race = factor(Race),
               Baseline.ECOG.Score = factor(Baseline.ECOG.Score),
               Met.Disease.Status = factor(Met.Disease.Status),
               Sample.age = factor(Sample.age),
               Received.platinum = factor(Received.platinum),
               Sample.collected.pre.platinum = factor(Sample.collected.pre.platinum),
               Lund = factor(Lund),
               Lund2 = factor(Lund2),
               Immune.phenotype = factor(Immune.phenotype),
               b8t = factor(b8t, levels = c("hi.lo", "lo.lo", "lo.hi", "hi.hi"))) %>% 
        select(-sample, -subject_id, -contains("_null"), -ENA.CHECKLIST, 
               -sizeFactor, -tag, -val, -binaryResponse, 
               -Best.Confirmed.Overall.Response, -Enrollment.IC, -TC.Level, -b_cell,
               -cd8_rose, -cd8_prat, -cd8_fehr, -ifng) %>% 
        relocate(Lund2, .before = "Lund")


# Univariate ---------------------------------------------------------------


# Logrank Test -------------------------------------------------------------

# Cutoff for consideration in multivariate analysis ("denoted 'sig'"): 
# p < 0.15

lr <- function(var) {
        name <- var
        var <- clin[[var]]
        coxph(Surv(os, censOS) ~ var, data = clin) %>% 
                glance() %>% 
                mutate(name = name) %>% 
                relocate(name)
}

vars <- colnames(clin)

vars <- vars[-which(vars %in% c("os", "censOS"))]

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
        filter(p.value.sc < 0.15) %>% 
        unite(term_2, name, stars, remove = F, sep = " ") %>% 
        select(-stars)

# Did not reach SS (p < 0.15):
# gender
# race
# intravesical bcg administered
# Tobacco Use History
# Tissue
# TCGA subtype


# Cox with Levels ---------------------------------------------------------

# Now look at coxph ratios for each level that reached ss

wl <- function(var) {
        feature <- var
        var <- clin[[var]]
        coxph(Surv(os, censOS) ~ var, data = clin) %>% 
                tidy(exponentiate = T) %>%  
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
                      "SE" = std.error,
                      "p-value" = p.value,
                      " " = stars) %>%
        select(-feature, -p.value.sc) %>% 
        relocate("p-value", .before =  " ") %>% 
        mutate(term_2 = str_replace_all(term_2, "\\.", " "),
               term_2 = str_replace(term_2,"IC Level", "PDL1 IC Level"),
               term_2 = str_replace(term_2, "FMOne mutation burden per MB", "Mutation Burden/MB"),
               term_2 = str_replace(term_2, "Sample age", "Sample Age"),
               term_2 = str_replace(term_2, "Received platinum", "Received Platinum"),
               term_2 = str_replace(term_2, "Sample collected pre platinum", "Sample Collected Pre-Platinum"),
               term_2 = str_replace(term_2, "Neoantigen burden per MB", "Neoantigen Burden/MB"),
               term_2 = str_replace(term_2, "Immune phenotype", "Immune Phenotype")) %>% 
        gt(groupname_col = "term_2",
           rowname_col = "term") %>% 
        fmt_number(columns = 2:3, drop_trailing_zeros = T) %>% 
        fmt_scientific(columns = 4, decimals = 2) %>% 
        cols_align("center") %>% 
        tab_style(style = cell_text(align = "right"), 
                  locations = cells_stub()) %>% 
        fmt_missing(columns = 1:6, missing_text = "") %>%
        tab_header(title = "Univariate Hazard Ratios") %>% 
        tab_footnote("Features not trending significant (P < 0.15) by logrank test: Sex, Race, Intravesical BCG, Tobacco Use History, TCGA Subtype, Tissue Site",
                     locations = cells_column_labels("p-value")) %>% 
        tab_footnote("p-value by logrank test (*** < 0.001, ** < 0.01, * < 0.05, # < 0.15)",
                     location = cells_row_groups(starts_with("PDL1"))) %>% 
        cols_width(vars(`p-value`) ~ px(130),
                   vars(HR, SE, " ") ~ px(60)) %>%
        tab_options(table.width = px(500)) %>% 
        gtsave("./figures/tables/risk_table_uni.png") 


# Multivariate -------------------------------------------------------

retidy <- function(df) {
        df %>% 
                mutate(feature = str_extract(term, colnames_search),
                       term = str_remove(term, colnames_search)) %>% 
                group_by(feature) %>% 
                select(-statistic)
}

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

mv <- coxph(Surv(os, censOS) ~ b8t + IC.Level + Neoantigen.burden.per.MB + 
                    Baseline.ECOG.Score + Met.Disease.Status + Received.platinum + 
                    Lund2, data = clin) %>% 
        tidy(exponentiate = T) %>%
        retidy() %>% 
        group_by(feature) %>% 
        group_split() %>% 
        lapply(make_reflevel_table) %>% 
        bind_rows() %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "(***)",
                                 p.value < 0.01 ~ "(**)",
                                 p.value < 0.05 ~ "(*)",
                                 T ~ "(NS)"))  %>% 
        filter(feature != "Received.platinum")

# Not significant:
# Received platinum


# Make Table ---------------------------------------------------------

mv %>% 
        dplyr::rename("HR" = estimate,
                      "SE" = std.error,
                      "p-value" = p.value,
                      " " = stars) %>%
        mutate(feature = str_replace_all(feature, "\\.", " "),
               feature = str_replace(feature, "FMOne mutation burden per MB", "Mutation Burden/MB"),
               feature = str_replace(feature, "Received platinum", "Received Platinum"),
               feature = str_replace(feature, "Sample collected pre platinum", "Sample Collected Pre-Platinum"),
               feature = str_replace(feature, "Neoantigen burden per MB", "Neoantigen Burden/MB"),
               feature = str_replace(feature, "b8t", "B8T"),
               term = str_replace(term, "hi.lo", "Hi/Lo"),
               term = str_replace(term, "lo.lo", "Lo/Lo"),
               term = str_replace(term, "lo.hi", "Lo/Hi"),
               term = str_replace(term, "hi.hi", "Hi/Hi")) %>% 
        gt(groupname_col = "feature",
           rowname_col = "term") %>% 
        fmt_number(columns = 2:4, drop_trailing_zeros = T) %>% 
        fmt_scientific(columns = 5, decimals = 2) %>% 
        cols_align("center") %>% 
        tab_style(style = cell_text(align = "right"), 
                  locations = cells_stub()) %>% 
        fmt_missing(columns = 1:6, missing_text = "") %>%
        tab_header(title = "Multivariate Hazard Ratios") %>% 
        tab_footnote("Feature not significant (P < 0.05) by logrank test: Received Platinum",
                     locations = cells_column_labels("p-value")) %>% 
        cols_width(vars(`p-value`) ~ px(130),
                   vars(HR, SE, " ") ~ px(60)) %>%
        tab_options(table.width = px(500)) %>% 
        gtsave("./figures/tables/risk_table_multi.png") 

