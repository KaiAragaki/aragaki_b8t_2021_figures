
# Description -------------------------------------------------------------

# Supplemental Figure 6


# Prepare Workspace --------------------------------------------------

library(tidyverse)
library(DESeq2)
library(survival)
library(survminer)
library(viridis)


# Combine datasets

riaz <- read_rds("../johnson_sex-diff_2020_figures/data/riaz/05_merge-colData.rds") %>% 
        filter(pre_or_on_treat == "Pre") %>% 
        select(Response, patient, Cohort, Dead.Alive...Dead...True., Time.to.Death...weeks., b_cell_spec:cd8_t_eff_null) %>% 
        mutate(Cohort = if_else(Cohort == "NIV3-NAIVE", F, T)) %>% 
        dplyr::rename(response = Response,
                      prior_ctla4 = Cohort,
                      dead = Dead.Alive...Dead...True.,
                      os = Time.to.Death...weeks.)

liu <- read_rds("../johnson_sex-diff_2020_figures/data/liu/04_merge-colData.Rds") %>% 
        filter(biopsyContext..1.Pre.Ipi..2.On.Ipi..3.Pre.PD1..4.On.PD1. == 3) %>% 
        select(sample, BR, OS, dead, priorCTLA4, b_cell_spec:cd8_t_eff_null) %>% 
        mutate(OS = OS/7,
               dead = as.logical(dead),
               priorCTLA4 = as.logical(priorCTLA4)) %>% 
        dplyr::rename(patient = sample,
                      response = BR,
                      prior_ctla4 = priorCTLA4,
                      os = OS)

col_data <- bind_rows(riaz, liu)


col_data <- col_data %>% 
        mutate(t_bin = if_else(cd8_t_eff > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               b_bin_spec = if_else(b_cell_spec > 0, "hi", "lo"),
               t_bin_spec = if_else(cd8_t_eff_spec > 0, "hi", "lo")) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        unite(b8t_spec, b_bin_spec, t_bin_spec, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               b8t_spec = factor(b8t_spec, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               response = factor(response, levels = c("CR", "PR", "SD", "PD", "MR", "NE")))



# Fig S6a: All Tumors vs Survival Across BCGS -----------------------------

png(filename = "./figures/fig_s6/fig_s6a.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ b_bin_spec, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.title = "B-cell Signature", 
                   legend.labs = c("Hi", "Lo"),
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 240),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()

# Fig S6b: All Tumors vs Survival Across T8GS -----------------------------

png(filename = "./figures/fig_s6/fig_s6b.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ t_bin_spec, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.title = "CD8+ T-cell Signature", 
                   legend.labs = c("Hi", "Lo"),
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 240),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Filter for just CTLA4 Naive ---------------------------------------------

naive <- col_data[col_data$prior_ctla4 == F, ]


# Fig S6c: Anti-CTLA4 Tumors vs Survival across CD8 T Signature -----------

png(filename = "./figures/fig_s6/fig_s6c.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ t_bin_spec, data = naive) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "CD8+ T-cell Signature\nCTLA4 Naive", 
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 240),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Filter for just Prior CTLA4 ----------------------------------------

pre <- col_data[col_data$prior_ctla4 == T, ]


# Fig S6c: Anti-CTLA4 Tumors vs Survival across CD8 T Signature -----------

png(filename = "./figures/fig_s6/fig_s6d.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ t_bin_spec, data = pre) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "CD8+ T-cell Signature\nPrior Anti-CTLA4", 
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 240),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()













png(filename = "./paper/figures/fig_4/riaz-liu/pre_t_surv_spec.png", width = 6, height = 4, units = "in", res = 288)
survfit(Surv(os, dead) ~ t_bin_spec, data = pre) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "CD8+ T-cell Signature", 
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 240),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()





