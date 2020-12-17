
# Description -------------------------------------------------------------

# Run GSVA


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(GSVA)
library(survival)
library(survminer)


# Read in Data ------------------------------------------------------------

exp <- read_rds("./data/merged/kim_sjodahl_choi_merged-exp.Rds")
clin <- read_tsv("./data/merged/kim_sjodahl_choi_merged-clin.tsv")
signatures <- read_rds("./data/signatures/signatures.Rds")

scores <- gsva(expr = exp, signatures, kcdf = "Gaussian", mx.diff = T)

plotting_dat <- t(scores) %>% 
        as_tibble(rownames = "geo_accession")

clin <- clin %>% 
        left_join(plotting_dat, by = "geo_accession") %>% 
        mutate(b_bin = if_else(b_cell > 0, "hi", "lo"),
               t_bin = if_else(cd8_rose > 0, "hi", "lo")) 

ggplot(plotting_dat, aes(b_cell)) + geom_density()
ggplot(plotting_dat, aes(cd8_rose)) + geom_density()
ggplot(clin, aes(x = b_cell, y = cd8_rose, color = dataset)) + geom_point()

survfit(Surv(clin$survival_month, clin$censor_os) ~ sex, data = clin) %>% 
        ggsurvplot(pval = T)

survfit(Surv(clin$survival_month, clin$censor_os) ~ t_bin, data = clin) %>% 
        ggsurvplot(pval = T)

survfit(Surv(clin$survival_month, clin$censor_os) ~ b_bin, data = clin) %>% 
        ggsurvplot(pval = T)

survfit(Surv(clin$survival_month, clin$censor_os) ~ t_bin + b_bin, data = clin) %>% 
        ggsurvplot(pval = T)

survfit(Surv(clin$survival_month, clin$censor_os) ~ t_bin + b_bin + sex, data = clin) %>% 
        ggsurvplot(pval = T)

survfit(Surv(clin$survival_month, clin$censor_os) ~ t_bin + sex, data = clin) %>% 
        ggsurvplot(pval = T)

survfit(Surv(clin$survival_month, clin$censor_os) ~ b_bin + sex, data = clin) %>% 
        ggsurvplot(pval = T)

survfit(Surv(clin$survival_month, clin$censor_os) ~ tumor, data = clin) %>% 
        ggsurvplot(pval = T)



female <- clin %>% 
        filter(sex == "F")

survfit(Surv(survival_month, censor_os) ~ t_bin + b_bin, data = female) %>% 
        ggsurvplot(pval = T)


male <- clin %>% 
        filter(sex == "M")

survfit(Surv(survival_month, censor_os) ~ t_bin + b_bin, data = male) %>% 
        ggsurvplot(pval = T)

