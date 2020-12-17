
# Description -------------------------------------------------------------

# Figure 6: The BCGS and the B8T signature segregates OS outcomes in
# patients treated with anti-PD-1 therapy for advanced melanoma in patients
# who have received prior CTLA-4


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(survival)
library(survminer)
library(viridis)


# Read in Data ------------------------------------------------------------


# Combine datasets
riaz <- read_rds("./data/riaz/05_merge-colData.rds") %>% 
        filter(pre_or_on_treat == "Pre") %>% 
        select(Response, 
               patient, 
               Cohort, 
               Dead.Alive...Dead...True., 
               Time.to.Death...weeks., 
               b_cell_spec:cd8_rose_null) %>% 
        mutate(Cohort = if_else(Cohort == "NIV3-NAIVE", F, T)) %>% 
        dplyr::rename(response = Response,
                      prior_ctla4 = Cohort,
                      dead = Dead.Alive...Dead...True.,
                      os = Time.to.Death...weeks.)

liu <- read_rds("./data/liu/04_merge-colData.Rds") %>% 
        filter(biopsyContext..1.Pre.Ipi..2.On.Ipi..3.Pre.PD1..4.On.PD1. == 3) %>% 
        select(sample, BR, OS, dead, priorCTLA4, b_cell_spec:cd8_rose_null) %>% 
        mutate(OS = OS/7,
               dead = as.logical(dead),
               priorCTLA4 = as.logical(priorCTLA4)) %>% 
        dplyr::rename(patient = sample,
                      response = BR,
                      prior_ctla4 = priorCTLA4,
                      os = OS)

col_data <- bind_rows(riaz, liu)

col_data <- col_data %>% 
        mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               b_bin_spec = if_else(b_cell_spec > 0, "hi", "lo"),
               t_bin_spec = if_else(cd8_rose_spec > 0, "hi", "lo")) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        unite(b8t_spec, b_bin_spec, t_bin_spec, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               b8t_spec = factor(b8t_spec, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               response = factor(response, levels = c("CR", "PR", "SD", "PD", "MR", "NE")))


# Fig 6a: All Patients B8T vs Survival ------------------------------------

png(filename = "./figures/fig_6/fig_6a.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ b8t_spec, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.title = "B8T Signature", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 228),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Filter Naive ------------------------------------------------------------

naive <- col_data[col_data$prior_ctla4 == F, ]


# Fig 6b: B-cell Signature vs Survival in Anti-CTLA4 Naive ----------------

png(filename = "./figures/fig_6/fig_6b.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ b_bin_spec, data = naive) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "B-cell Signature",
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 228),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig 6c: B8T vs Survival in Anti-CTLA4 Naive -----------------------------

png(filename = "./figures/fig_6/fig_6c.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ b8t_spec, data = naive) %>% 
        ggsurvplot(pval = T,
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.title = "B8T Signature", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 228),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()



# Filter Prior ------------------------------------------------------------

pre <- col_data[col_data$prior_ctla4 == T, ]


# Fig 6d: B-cell Signature vs Survival in Anti-CTLA4 Prior ----------------

png(filename = "./figures/fig_6/fig_6d.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ b_bin_spec, data = pre) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "B-cell Signature",
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 228),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig 6e: B-cell Signature vs Survival in Anti-CTLA4 Prior ----------------

fig_6e <- filter(pre, b8t_spec != "hi_lo")

png(filename = "./figures/fig_6/fig_6e.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, dead) ~ b8t_spec, data = fig_6e) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.25,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "B8T Signature", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo"),
                   xlab = "Overall Survival (Months)",
                   legend = "right",
                   xscale = 4,
                   break.x.by = 40,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 228),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()
