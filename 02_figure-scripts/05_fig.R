
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
               b8t_bin = if_else(b_cell>0 & cd8_rose>0,"hi_hi","all_other")) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               patient.gender = factor(patient.gender, levels = c("male", "female")),
               b8t_bin = factor(b8t_bin, levels=c("hi_hi","all_other")))


# Fig 5a: B8T vs TCGA BLCA Survival ---------------------------------------

png(filename = "./figures/fig_5/fig_5a.png", width = 5.5, height = 5.7, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1),
                   risk.table = T,
                   title = "A.",
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   #subtitle="A    ALL MIBC Tumors",
                   legend = "right",
                   xscale = 30,
                   break.x.by = 300,
                   xlim = c(0, 1800),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


<<<<<<< HEAD

# Fig 5a2: Hi/Hi vs Others -------------------------------------------
=======
# Fig 5b: Hi/Hi vs Others -------------------------------------------
>>>>>>> 2b480769ac0dace750ae34f138d4835690a4e216

png(filename = "./figures/fig_5/fig_5b.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_bin, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1),
                   risk.table = T,
                   title = "B.",
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T", 
                   legend.labs = c("Hi/Hi", "Other"),
                   xlab = "Overall Survival (Months)", 
                   #subtitle="A2    ALL MIBC Tumors",
                   legend ="right",
                   xscale = 30,
                   break.x.by = 300,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   xlim = c(0, 1800),
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


# Fig 5c: B8T vs TCGA BLCA Survival in Females -----------------------------

fig_5c <- col_data %>% 
        filter(patient.gender == "female")

png(filename = "./figures/fig_5/fig_5c.png", width = 5.5, height = 5.7, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t, data = fig_5c) %>% 
        ggsurvplot(pval = T,
                   palette = viridis(4, end = 0.95, direction = -1),
                   risk.table = T,
                   title = "C.",
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T\n(Females)", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   #subtitle="B     Female MIBC",
                   legend = "right",
                   xscale = 30,
                   break.x.by = 300,
                   xlim = c(0, 1800),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


<<<<<<< HEAD

# Fig 5b2: Hi/Hi vs Others - Females ---------------------------------
=======
# Fig 5d: Hi/Hi vs Others - Females ---------------------------------
>>>>>>> 2b480769ac0dace750ae34f138d4835690a4e216

png(filename = "./figures/fig_5/fig_5d.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_bin, data = fig_5c) %>% 
        ggsurvplot(pval = T,
                   palette = viridis(2, end = 0.95, direction = -1),
                   risk.table = T,
                   title = "D.",
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T\n(Females)", 
<<<<<<< HEAD
                   legend.labs = c("Hi/Hi", "All other"),
                   xlab = "Overall Survival (Months)",
                   #subtitle="B2     Female MIBC",
=======
                   legend.labs = c("Hi/Hi", "Other"),
                   xlab = "Overall Survival (Months)", 
>>>>>>> 2b480769ac0dace750ae34f138d4835690a4e216
                   legend = "right",
                   xscale = 30,
                   break.x.by = 300,
                   xlim = c(0, 1800),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


<<<<<<< HEAD

# Fig 5c: B8T vs TCGA BLCA Survival in Males -------------------------------
=======
# Fig e: B8T vs TCGA BLCA Survival in Males -------------------------------
>>>>>>> 2b480769ac0dace750ae34f138d4835690a4e216

fig_5e <- col_data %>% 
        filter(patient.gender == "male")

png(filename = "./figures/fig_5/fig_5e.png", width = 5.5, height = 5.7, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t, data = fig_5e) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1),
                   risk.table = T,
                   title = "E.",
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T\n(Males)", 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   #subtitle="C     Male MIBC",
                   legend = "right",
                   xscale = 30,
                   break.x.by = 300,
                   xlim = c(0, 1800),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


# Fig 5f: Hi/Hi vs Others - Males -----------------------------------

png(filename = "./figures/fig_5/fig_5f.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_bin, data = fig_5e) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1),
                   risk.table = T,
                   title = "F.",
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   legend.title = "B8T\n(Males)", 
<<<<<<< HEAD
                   legend.labs = c("Hi/Hi", "All other"),
                   xlab = "Overall Survival (Months)",
                   #subtitle="C2     Male MIBC",
=======
                   legend.labs = c("Hi/Hi", "Other"),
                   xlab = "Overall Survival (Months)", 
>>>>>>> 2b480769ac0dace750ae34f138d4835690a4e216
                   legend = "right",
                   xscale = 30,
                   break.x.by = 300,
                   xlim = c(0, 1800),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()

