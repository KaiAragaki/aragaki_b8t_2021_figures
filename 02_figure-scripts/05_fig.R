
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
        mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"), 
               b8t_bin = if_else(b_cell > 0 & cd8_rose > 0,"hi_hi", "all_other")) %>%
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               patient.gender = factor(patient.gender, levels = c("male", "female")),
               b8t_bin = factor(b8t_bin, levels=c("hi_hi","all_other")))


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


# Fig 5a2: Hi/Hi vs Others -------------------------------------------

png(filename = "./figures/fig_5/fig_5a2.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_bin, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1),
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T", 
                   legend.labs = c("Hi/Hi", "All other"),
                   xlab = "Overall Survival (Months)", 
                   legend ="right",
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
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t, data = fig_5b) %>% 
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


# Fig 5b2: Hi/Hi vs Others - Females ---------------------------------

png(filename = "./figures/fig_5/fig_5b2.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_bin, data = fig_5b) %>% 
        ggsurvplot(pval = T,
                   palette = viridis(2, end = 0.95, direction = -1),
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T\n(Females)", 
                   legend.labs = c("Hi/Hi", "All other"),
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



# Fig 5c2: Hi/Hi vs Others - Males -----------------------------------


png(filename = "./figures/fig_5/fig_5c2.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_bin, data = fig_5c) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1),
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T\n(Males)", 
                   legend.labs = c("Hi/Hi", "All other"),
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