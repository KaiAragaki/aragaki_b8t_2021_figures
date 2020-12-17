
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
               t.bin = if_else(cd8_t_eff > 0, "hi", "lo")) %>% 
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
               Immune.phenotype = factor(Immune.phenotype)) %>% 
        relocate(Lund2, .before = "Lund")


# Univariate ---------------------------------------------------------


# Anovas -------------------------------------------------------------



# Cutoff for consideration in multivariate analysis ("denoted 'sig'"): 
# p < 0.15

a <- coxph(Surv(os, censOS)~ gender, data = clin) %>% 
        anova() %>%  
        tidy()
b <- coxph(Surv(os, censOS)~IC.Level, data = clin) %>% 
        anova() %>% 
        tidy()
c <- coxph(Surv(os, censOS)~FMOne.mutation.burden.per.MB, data = clin) %>% 
        anova() %>% 
        tidy()
d <- coxph(Surv(os, censOS)~Race, data = clin) %>% 
        anova() %>% 
        tidy()
e <- coxph(Surv(os, censOS)~Intravesical.BCG.administered, data = clin) %>% 
        anova() %>% 
        tidy()
f <- coxph(Surv(os, censOS)~Baseline.ECOG.Score, data = clin) %>% 
        anova() %>% 
        tidy()
g <- coxph(Surv(os, censOS)~Tobacco.Use.History, data = clin) %>% 
        anova() %>% 
        tidy()
h <- coxph(Surv(os, censOS)~Met.Disease.Status, data = clin) %>% 
        anova() %>% 
        tidy()
i <- coxph(Surv(os, censOS)~Sample.age, data = clin) %>% 
        anova() %>% 
        tidy()
j <- coxph(Surv(os, censOS)~Received.platinum, data = clin) %>% 
        anova() %>% 
        tidy()
k <- coxph(Surv(os, censOS)~Sample.collected.pre.platinum, data = clin) %>% 
        anova() %>% 
        tidy()
l <- coxph(Surv(os, censOS)~Neoantigen.burden.per.MB, data = clin) %>% 
        anova() %>% 
        tidy()
m <- coxph(Surv(os, censOS)~Lund, data = clin) %>% 
        anova() %>% 
        tidy()
n <- coxph(Surv(os, censOS)~Lund2, data = clin) %>% 
        anova() %>% 
        tidy()
o <- coxph(Surv(os, censOS)~TCGA.Subtype, data = clin) %>% 
        anova() %>% 
        tidy()
p <- coxph(Surv(os, censOS)~Immune.phenotype, data = clin) %>% 
        anova() %>% 
        tidy()

signif_uni <- rbind(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) %>% 
        filter(term != "NULL") %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "(***)",
                                 p.value < 0.01 ~ "(**)",
                                 p.value < 0.05 ~ "(*)",
                                 p.value < 0.15 ~ "(#)",
                                 T ~ "(NS)")) %>% 
        filter(p.value < 0.15) %>% 
        unite(term_2, term, stars, sep = " ", remove = F)

# Did not reach SS (p < 0.15):
# gender
# race
# intravesical bcg administered
# Tobacco Use History
# TCGA subtype




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



a.1 <- coxph(Surv(os, censOS)~IC.Level, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
b.1 <- coxph(Surv(os, censOS)~FMOne.mutation.burden.per.MB, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
c.1 <- coxph(Surv(os, censOS)~Baseline.ECOG.Score, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
d.1 <- coxph(Surv(os, censOS)~Met.Disease.Status, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
e.1 <- coxph(Surv(os, censOS)~Sample.age, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
f.1 <- coxph(Surv(os, censOS)~Received.platinum, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
g.1 <- coxph(Surv(os, censOS)~Sample.collected.pre.platinum, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
h.1 <- coxph(Surv(os, censOS)~Neoantigen.burden.per.MB, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
i.1 <- coxph(Surv(os, censOS)~Lund, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
j.1 <- coxph(Surv(os, censOS)~Lund2, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()
k.1 <- coxph(Surv(os, censOS)~Immune.phenotype, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        retidy()


# Make Table ---------------------------------------------------------

bind_rows(a.1, b.1, c.1, d.1, e.1, f.1, g.1, h.1, i.1, j.1, k.1) %>% 
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
        fmt_number(columns = 2:3) %>% 
        fmt_scientific(columns = 4, decimals = 2) %>% 
        cols_align("center") %>% 
        tab_style(style = cell_text(align = "right"), 
                  locations = cells_stub()) %>% 
        fmt_missing(columns = 1:6, missing_text = "") %>%
        tab_header(title = "Univariate Hazard Ratios") %>% 
        tab_footnote("Features not trending significant (P < 0.15) by ANOVA: Gender, Race, Intravesical BCG, Tobacco Use History, TCGA Subtype",
                     locations = cells_column_labels("p-value")) %>% 
        tab_footnote("p-value by ANOVA (*** < 0.001, ** < 0.01, * < 0.05, # < 0.15)",
                     location = cells_row_groups(starts_with("PDL1"))) %>% 
        cols_width(vars(`p-value`) ~ px(130),
                   vars(HR, SE, " ") ~ px(60)) %>%
        tab_options(table.width = px(500)) %>% 
        gtsave("./figures/tables/risk_table_uni.png") 




# Multivariate -------------------------------------------------------


# Anovas -------------------------------------------------------------


a <- coxph(Surv(os, censOS)~b8t + IC.Level, data = clin) %>% 
        anova() %>% 
        tidy()
b <- coxph(Surv(os, censOS)~b8t + FMOne.mutation.burden.per.MB, data = clin) %>% 
        anova() %>% 
        tidy()
c <- coxph(Surv(os, censOS)~b8t + Baseline.ECOG.Score, data = clin) %>% 
        anova() %>% 
        tidy()
d <- coxph(Surv(os, censOS)~b8t + Met.Disease.Status, data = clin) %>% 
        anova() %>% 
        tidy()
e <- coxph(Surv(os, censOS)~b8t + Sample.age, data = clin) %>% 
        anova() %>% 
        tidy()
f <- coxph(Surv(os, censOS)~b8t + Received.platinum, data = clin) %>% 
        anova() %>% 
        tidy()
g <- coxph(Surv(os, censOS)~b8t + Sample.collected.pre.platinum, data = clin) %>% 
        anova() %>% 
        tidy()
h <- coxph(Surv(os, censOS)~b8t + Neoantigen.burden.per.MB, data = clin) %>% 
        anova() %>% 
        tidy()
i <- coxph(Surv(os, censOS)~b8t + Lund, data = clin) %>% 
        anova() %>% 
        tidy()
j <- coxph(Surv(os, censOS)~b8t + Lund2, data = clin) %>% 
        anova() %>% 
        tidy()
k <- coxph(Surv(os, censOS)~b8t + Immune.phenotype, data = clin) %>% 
        anova() %>% 
        tidy()

signif_multi <- rbind(a, b, c, d, e, f, g, h, i, j, k) %>% 
        filter(term != "NULL") %>%
        filter(term != "b8t") %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "(***)",
                                 p.value < 0.01 ~ "(**)",
                                 p.value < 0.05 ~ "(*)",
                                 T ~ "(NS)")) %>% 
        filter(p.value < 0.05) %>% 
        unite(term_2, term, stars, sep = " ", remove = F)

# Not significant:
# Immune phenotype
# Sample age
# IC Level


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


a.1 <- coxph(Surv(os, censOS)~b8t + FMOne.mutation.burden.per.MB, data = clin) %>% 
        tidy(exponentiate = T) %>%
        filter(str_detect(term, "^b8t", negate = T)) %>% 
        retidy()
b.1 <- coxph(Surv(os, censOS)~b8t + Baseline.ECOG.Score, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        filter(str_detect(term, "^b8t", negate = T)) %>% 
        retidy()
c.1 <- coxph(Surv(os, censOS)~b8t + Met.Disease.Status, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        filter(str_detect(term, "^b8t", negate = T)) %>% 
        retidy()
d.1 <- coxph(Surv(os, censOS)~b8t + Received.platinum, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        filter(str_detect(term, "^b8t", negate = T)) %>% 
        retidy()
e.1 <- coxph(Surv(os, censOS)~b8t + Sample.collected.pre.platinum, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        filter(str_detect(term, "^b8t", negate = T)) %>% 
        retidy()
f.1 <- coxph(Surv(os, censOS)~b8t + Neoantigen.burden.per.MB, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        filter(str_detect(term, "^b8t", negate = T)) %>% 
        retidy()
g.1 <- coxph(Surv(os, censOS)~b8t + Lund, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        filter(str_detect(term, "^b8t", negate = T)) %>% 
        retidy()
h.1 <- coxph(Surv(os, censOS)~b8t + Lund2, data = clin) %>% 
        tidy(exponentiate = T) %>% 
        filter(str_detect(term, "^b8t", negate = T)) %>% 
        retidy()


# Make Table ---------------------------------------------------------

bind_rows(a.1, b.1, c.1, d.1, e.1, f.1, g.1, h.1) %>% 
        mutate(stars = case_when(p.value < 0.001 ~ "***",
                                 p.value < 0.01 ~ "**",
                                 p.value < 0.05 ~ "*",
                                 is.na(p.value) ~ NA_character_,
                                 T ~ "NS"))%>% 
        dplyr::rename("HR" = estimate,
                      "SE" = std.error,
                      "p-value" = p.value,
                      " " = stars) %>%
        left_join(signif_multi, by = c("feature" = "term")) %>%
        select(-(logLik:stars), -feature) %>% 
        mutate(term_2 = str_replace_all(term_2, "\\.", " "),
               term_2 = str_replace(term_2, "FMOne mutation burden per MB", "Mutation Burden/MB"),
               term_2 = str_replace(term_2, "Received platinum", "Received Platinum"),
               term_2 = str_replace(term_2, "Sample collected pre platinum", "Sample Collected Pre-Platinum"),
               term_2 = str_replace(term_2, "Neoantigen burden per MB", "Neoantigen Burden/MB")) %>% 
        gt(groupname_col = "term_2",
           rowname_col = "term") %>% 
        fmt_number(columns = 2:3) %>% 
        fmt_scientific(columns = 4, decimals = 2) %>% 
        cols_align("center") %>% 
        tab_style(style = cell_text(align = "right"), 
                  locations = cells_stub()) %>% 
        fmt_missing(columns = 1:6, missing_text = "") %>%
        tab_header(title = "Multivariate Hazard Ratios") %>% 
        tab_footnote("Features not significant (P < 0.05) by ANOVA: Immune Phenotype, Sample Age, IC Level",
                     locations = cells_column_labels("p-value")) %>% 
        tab_footnote("p-value by ANOVA (*** < 0.001, ** < 0.01, * < 0.05)",
                     location = cells_row_groups(starts_with("Mutation"))) %>% 
        cols_width(vars(`p-value`) ~ px(130),
                   vars(HR, SE, " ") ~ px(60)) %>%
        tab_options(table.width = px(500)) %>% 
        gtsave("./figures/tables/risk_table_multi.png") 


# Adjusted Curves ----------------------------------------------------

b8t <- clin %>% 
        unite(b8t, b.bin_, t.bin_) %>% 
        as.data.frame()

ic <- coxph(Surv(os_, censOS_)~strata(b8t) + IC.Level_, b8t)
ecog <- coxph(Surv(os_, censOS_)~strata(b8t) + Baseline.ECOG.Score_, b8t)
mut <- coxph(Surv(os_, censOS_)~strata(b8t) + FMOne.mutation.burden.per.MB_, b8t)
met <- coxph(Surv(os_, censOS_)~strata(b8t) + Met.Disease.Status_, b8t)
samp_age <- coxph(Surv(os_, censOS_)~strata(b8t) + Sample.age_, b8t)
plat <- coxph(Surv(os_, censOS_)~strata(b8t) + Received.platinum_, b8t)
is_preplat <- coxph(Surv(os_, censOS_)~strata(b8t) + Sample.collected.pre.platinum_, b8t)
mut_neo <- coxph(Surv(os_, censOS_)~strata(b8t)+ Neoantigen.burden.per.MB_, b8t)
lund <- coxph(Surv(os_, censOS_)~strata(b8t) + Lund_, b8t)
lund2 <- coxph(Surv(os_, censOS_)~strata(b8t) + Lund2_, b8t)
immune <- coxph(Surv(os_, censOS_)~strata(b8t) + Immune.phenotype_, b8t)


ggadjustedcurves(og)
ggadjustedcurves(ic)
ggadjustedcurves(ec)
ggadjustedcurves(mut)
ggadjustedcurves(met)
ggadjustedcurves(samp_age)
ggadjustedcurves(plat)
ggadjustedcurves(is_preplat)
ggadjustedcurves(mut_neo)
ggadjustedcurves(lund)
ggadjustedcurves(lund2)
ggadjustedcurves(immune)


# ----- BUILDING A COXPH MODEL -----
# Fwd selection, forcing B8T as first --------------------------------

aaa <- coxph(Surv(os_, censOS_)~b8t_, clin)
a <- coxph(Surv(os_, censOS_)~b8t_ + IC.Level_, clin)
b <- coxph(Surv(os_, censOS_)~b8t_ + FMOne.mutation.burden.per.MB_, clin)
c <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_, clin)
d <- coxph(Surv(os_, censOS_)~b8t_ + Met.Disease.Status_, clin)
e <- coxph(Surv(os_, censOS_)~b8t_ + Sample.age_, clin)
f <- coxph(Surv(os_, censOS_)~b8t_ + Received.platinum_, clin)
g <- coxph(Surv(os_, censOS_)~b8t_ + Sample.collected.pre.platinum_, clin)
h <- coxph(Surv(os_, censOS_)~b8t_ + Neoantigen.burden.per.MB_, clin)
i <- coxph(Surv(os_, censOS_)~b8t_ + Lund_, clin)
j <- coxph(Surv(os_, censOS_)~b8t_ + Lund2_, clin)
k <- coxph(Surv(os_, censOS_)~b8t_ + Immune.phenotype_, clin)
l <- coxph(Surv(os_, censOS_)~b8t_ + gender_, clin)

anova(a)
anova(b) # 1.04e-4
anova(c) # 1.37e-6 <- lowest - ECOG
anova(d) # 2.21e-6
anova(e) 
anova(f) # 0.04594
anova(g) # 0.0252
anova(h) # 8.91e-5
anova(i) # 0.03299
anova(j) # 0.01048
anova(k)
anova(l)


a.1 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + FMOne.mutation.burden.per.MB_, clin)
b.1 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Met.Disease.Status_, clin)
c.1 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Received.platinum_, clin)
d.1 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Sample.collected.pre.platinum_, clin)
e.1 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_, clin)
f.1 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Lund_, clin)
g.1 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Lund2_, clin)

anova(a.1) # 3.45e-6
anova(b.1) # 6.06e-5
anova(c.1) # 1.48e-3
anova(d.1)
anova(e.1) # 1.14e-6
anova(f.1) # 1.04e-2
anova(g.1) # 2.53e-3


a.2 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + FMOne.mutation.burden.per.MB_, clin)
b.2 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Met.Disease.Status_, clin)
c.2 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Received.platinum_, clin)
d.2 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund_, clin)
e.2 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund2_, clin)

anova(a.2)
anova(b.2) # 7.14e-3
anova(c.2)
anova(d.2)
anova(e.2) # 0.04


a.3 <- coxph(Surv(os_, censOS_)~b8t_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Met.Disease.Status_ + Lund2_, clin, )

anova(a.3)   # 0.004

# Fwd selection, withholding B8T until end ---------------------------

a <- coxph(Surv(os_, censOS_)~IC.Level_, clin)
b <- coxph(Surv(os_, censOS_)~FMOne.mutation.burden.per.MB_, clin)
c <- coxph(Surv(os_, censOS_)~Baseline.ECOG.Score_, clin)
d <- coxph(Surv(os_, censOS_)~Met.Disease.Status_, clin)
e <- coxph(Surv(os_, censOS_)~Sample.age_, clin)
f <- coxph(Surv(os_, censOS_)~Received.platinum_, clin)
g <- coxph(Surv(os_, censOS_)~Sample.collected.pre.platinum_, clin)
h <- coxph(Surv(os_, censOS_)~Neoantigen.burden.per.MB_, clin)
i <- coxph(Surv(os_, censOS_)~Lund_, clin)
j <- coxph(Surv(os_, censOS_)~Lund2_, clin)
k <- coxph(Surv(os_, censOS_)~Immune.phenotype_, clin)
l <- coxph(Surv(os_, censOS_)~gender_, clin)

anova(a) # 1.87e-3
anova(b) # 1.34e-5
anova(c) # 1.54e-6
anova(d) # 1.19e-6
anova(e) # 0.10
anova(f) # 3.90e-2
anova(g) # 6.30e-2
anova(h) # 4.54e-6
anova(i) # 6.84e-2
anova(j) # 1.98e-2
anova(k) # 8.06e-2
anova(l)


a.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + IC.Level_, clin)
b.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + FMOne.mutation.burden.per.MB_, clin)
c.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_, clin)
d.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Sample.age_, clin)
e.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Received.platinum_, clin)
f.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Sample.collected.pre.platinum_, clin)
g.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Neoantigen.burden.per.MB_, clin)
h.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Lund_, clin)
i.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Lund2_, clin)
j.1 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Immune.phenotype_, clin)

anova(a.1) # 1.04e-2
anova(b.1) # 2.40e-3
anova(c.1) # 6.53e-5
anova(d.1) 
anova(e.1) 
anova(f.1)
anova(g.1) # 3.17e-4
anova(h.1) # 1.35e-2
anova(i.1) # 4.88e-3
anova(j.1)


a.2 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + IC.Level_, clin)
b.2 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + FMOne.mutation.burden.per.MB_, clin)
c.2 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_, clin)
d.2 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Lund_, clin)
e.2 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Lund2_, clin)

anova(a.2) # 7.24e-3
anova(b.2) # 7.83e-5
anova(c.2) # 1.85e-5
anova(d.2) # 9.57e-3
anova(e.2) # 3.14e-3


a.3 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + IC.Level_, clin)
b.3 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + FMOne.mutation.burden.per.MB_, clin)
c.3 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund_, clin)
d.3 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund2_, clin)

anova(a.3) # 3.90e-2
anova(b.3)
anova(c.3) # 3.07e-2
anova(d.3) # 8.75e-3

a.4 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund2_ + IC.Level_, clin)
b.4 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund2_ + Lund_, clin)

anova(a.4) # 3.67e-3
anova(b.4)

a.5 <- coxph(Surv(os_, censOS_)~Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund2_ + IC.Level_ + b8t_, clin)

anova(a.5)


# Fwd selection, B8T as predictor, not forcing first -----------------

a <- coxph(Surv(os_, censOS_)~IC.Level_, clin)
b <- coxph(Surv(os_, censOS_)~FMOne.mutation.burden.per.MB_, clin)
c <- coxph(Surv(os_, censOS_)~Baseline.ECOG.Score_, clin)
d <- coxph(Surv(os_, censOS_)~Met.Disease.Status_, clin)
e <- coxph(Surv(os_, censOS_)~Sample.age_, clin)
g <- coxph(Surv(os_, censOS_)~Sample.collected.pre.platinum_, clin)
h <- coxph(Surv(os_, censOS_)~Neoantigen.burden.per.MB_, clin)
i <- coxph(Surv(os_, censOS_)~Lund_, clin)
j <- coxph(Surv(os_, censOS_)~Lund2_, clin)
k <- coxph(Surv(os_, censOS_)~Immune.phenotype_, clin)
l <- coxph(Surv(os_, censOS_)~gender_, clin)
m <- coxph(Surv(os_, censOS_)~b8t_, clin)

anova(a) # 1.87e-3
anova(b) # 1.34e-5
anova(c) # 1.54e-6
anova(d) # 1.19e-6
anova(e) # 0.10
anova(f) # 0.04
anova(g) # 0.06
anova(h) # 4.54e-6
anova(i) # 0.07
anova(j) # 0.02
anova(k) # 0.08
anova(l) # 0.18
anova(m) # 3.62e-5


a.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + IC.Level_, clin)
b.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + FMOne.mutation.burden.per.MB_, clin)
c.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_, clin)
d.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Sample.age_, clin)
e.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Sample.collected.pre.platinum_, clin)
f.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Neoantigen.burden.per.MB_, clin)
g.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Lund_, clin)
h.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Lund2_, clin)
i.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Immune.phenotype_, clin)
j.1 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + b8t_, clin)

# Keep P < 0.15 from univariate setting
anova(a.1) # 0.01
anova(b.1) # 2.40e-3
anova(c.1) # 6.53e-5
anova(d.1) 
anova(e.1) 
anova(f.1) # 3.17e-4
anova(g.1) # 1.35e-2
anova(h.1) # 4.88e-3
anova(i.1) 
anova(j.1) # 9.53e-5

a.2 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + IC.Level_, clin)
b.2 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + FMOne.mutation.burden.per.MB_, clin)
c.2 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_, clin)
d.2 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Lund_, clin)
e.2 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Lund2_, clin)
f.2 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + b8t_, clin)

anova(a.2) # 7.24e-3
anova(b.2) # 7.83e-5
anova(c.2) # 1.85e-5
anova(d.2) # 9.57e-3
anova(e.2) # 3.14e-3
anova(f.2) # 7.83e-5

a.3 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + IC.Level_, clin)
b.3 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + FMOne.mutation.burden.per.MB_, clin)
c.3 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund_, clin)
d.3 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund2_, clin)
e.3 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + b8t_, clin)

anova(a.3) # 3.90e-2 
anova(b.3)
anova(c.3) # 3.07e-2
anova(d.3) # 8.75e-3
anova(e.3) # 9.67e-5

a.4 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + b8t_ + IC.Level_, clin)
b.4 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + b8t_ + Lund_, clin)
c.4 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + b8t_ + Lund2_, clin)

anova(a.4) 
anova(b.4) # 1.21e-2
anova(c.4) # 4.19e-3

a.6 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + b8t_ + Lund2_ + Lund_, clin)

anova(a.6)


# Summary ------------------------------------------------------------

# Best model if B8T is the first predictor:

a.3 <- coxph(Surv(os_, censOS_)~ b8t_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Met.Disease.Status_ + Lund2_, clin)

# Best model if B8T is withheld until end:

a.5 <- coxph(Surv(os_, censOS_)~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + Lund2_ + IC.Level_ + b8t_, clin)

# Best model if B8T is used as just another predictor:

a.6 <- coxph(Surv(os_, censOS_) ~ Met.Disease.Status_ + Baseline.ECOG.Score_ + Neoantigen.burden.per.MB_ + b8t_ + Lund2_, clin)

# All models converge with identical predictors

