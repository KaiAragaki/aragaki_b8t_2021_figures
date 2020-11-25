
# Description -------------------------------------------------------------

# Figure 5: The B8T signature segregates OS outcomes in non-checkpoint treated
# tumors for women, but not men.


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(survival)
library(survminer)
library(viridis)

blca <- read_rds("./data/tcga-blca/A_02_normalized-counts.Rds")
sig_scores <- read_rds("./data/tcga-blca/C_03_gsva-scores.Rds") %>%
        t()

if(all(colnames(blca) == rownames(sig_scores))){ # Make sure everything is arranged properly
        colData(blca) <- cbind(colData(blca), sig_scores)
        print(T)
}

col_data <- colData(blca) %>% 
        as_tibble() %>% 
        mutate(t_bin = if_else(cd8_t_eff > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo")) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               patient.gender = factor(patient.gender, levels = c("male", "female")))


# Fig 5a: B8T vs TCGA BLCA Survival ---------------------------------------

png(filename = "./figures/fig_5/fig_5a.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1),
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 30,
                   break.x.by = 300,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 1800),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig 5b: B8T vs TCGA BLCA Survival in Females -----------------------------

fig_5b <- col_data %>% 
        filter(patient.gender == "female")

png(filename = "./figures/fig_5/fig_5b.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t, data = mysample) %>% 
        ggsurvplot(pval = T,
                   palette = viridis(4, end = 0.95, direction = -1),
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T\n(Females)", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 30,
                   break.x.by = 300,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 1800),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()

# Fig 5c: B8T vs TCGA BLCA Survival in Males -------------------------------

fig_5c <- col_data %>% 
        filter(patient.gender == "male")

png(filename = "./figures/fig_5/fig_5c.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t, data = fig_5c) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1),
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T\n(Males)", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   legend = "right",
                   xscale = 30,
                   break.x.by = 300,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 1800),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()



# Read in Melanoma --------------------------------------------------------

col_data <- read_rds("./data/tcga-skcm/C_03_merge-colData.rds") %>% 
        mutate(t_bin = if_else(cd8_t_eff > 0, "hi", "lo"),
               t_bin_spec = if_else(cd8_t_eff_spec > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               b_bin_spec = if_else(b_cell > 0, "hi", "lo"),
               patient.gender = factor(patient.gender, c("male", "female"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        unite(b8t_spec, b_bin_spec, t_bin_spec, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               b8t_spec = factor(b8t_spec, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")))
lym <- filter(col_data, patient.tumor_tissue_site == "regional lymph node")


# Fig 5d: Female Lymph Node Metastases: B8T vs Survival -------------------

fig_5d <- filter(lym, patient.gender == "female")

png(filename = "./figures/fig_5/fig_5d.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_spec, data = fig_5d) %>% 
        ggsurvplot(pval = T,
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"), 
                   legend.title = "B8T\n(Sex = F)", 
                   xlab = "Overall Survival (Years)", 
                   legend = "right",
                   xscale = 365,
                   break.x.by = 730,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 6000),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig 5e: Male Lymph Node Metastases: B8T vs Survival ---------------------

fig_5e <- filter(lym, patient.gender == "male")

png(filename = "./figures/fig_5/fig_5e.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_spec, data = fig_5e) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"), 
                   legend.title = "B8T\n(Sex = M)", 
                   xlab = "Overall Survival (Years)", 
                   legend = "right",
                   xscale = 365,
                   break.x.by = 730,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 6000),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()

