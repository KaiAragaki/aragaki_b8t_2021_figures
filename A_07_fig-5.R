
# Description -------------------------------------------------------------

# Figure 5: The B8T signature segregates OS outcomes in non-checkpoint treated
# tumors for women, but not men.


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(survival)
library(survminer)
library(viridis)



blca <- read_rds("H:/Rusty_B8T/data/TCGA_BLCA/A_02_normalized-counts.Rds")

sig_scores <- read_rds("H:/Rusty_B8T/data/TCGA_BLCA/C_03_gsva-scores.Rds") %>%
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

dim(col_data)
View(col_data)

data_b8t=col_data%>%
          select(id,path,shortID,mRNA.cluster,patient.gender,death_days,followUp_days,new_death,death_event,b8t, t_bin,b_bin, b_cell,cd8_t_eff)%>%
  as_tibble()
          
View(data_b8t)
write_csv(data_b8t, "H:/Rusty_B8T/data/TCGA_BLCA/data_b8t.csv")


# Fig 5a: B8T vs TCGA BLCA Survival ---------------------------------------

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_5a.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t, data = col_data) %>% 
          ggsurvplot(pval = T, 
             palette = viridis(4, end = 0.95, direction = -1),
             risk.table = T,
             risk.table.height = 0.3,
             risk.table.title = "No. at risk",
             legend.title = "B-cell/Cd8+ T-cell", 
             legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
             xlab = "Overall Survival (Months)", 
             subtitle="A    ALL MIBC Tumors",
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
summary(ggsurv)
View(ggsurv$table)

# Fig 5b: B8T vs TCGA BLCA Survival in Females -----------------------------

fig_5b <- col_data %>% 
  filter(patient.gender == "female")

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_5b.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurvb <- survfit(Surv(new_death, death_event) ~ b8t, data = fig_5b) %>% 
  ggsurvplot(pval = T,
             palette = viridis(4, end = 0.95, direction = -1),
             risk.table = T,
             risk.table.height = 0.3,
             risk.table.title = "No. at risk",
             legend.title = "B-cell/Cd8+ T-cell\n(Females)", 
             legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
             xlab = "Overall Survival (Months)", 
             legend = "right",
             subtitle="B     Female MIBC",
             xscale = 30,
             break.x.by = 300,
             font.xtickslab = 15, 
             font.ytickslab = 15, 
             font.legend = 12,
             xlim = c(0, 1800),
             pval.coord = c(0, 0.05))
ggsurvb$table <- ggsurvb$table +
  ylab(NULL) + 
  xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsurvb
dev.off()

# Fig 5c: B8T vs TCGA BLCA Survival in Males -------------------------------

fig_5c <- col_data %>% 
  filter(patient.gender == "male")

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_5c.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurvc <- survfit(Surv(new_death, death_event) ~ b8t, data = fig_5c) %>% 
  ggsurvplot(pval = T, 
             palette = viridis(4, end = 0.95, direction = -1),
             risk.table = T,
             risk.table.height = 0.3,
             risk.table.title = "No. at risk",
             legend.title = "B-cell/Cd8+ T-cell\n(Males)", 
             legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
             xlab = "Overall Survival (Months)", 
             legend = "right",
             subtitle="C     Male MIBC",
             xscale = 30,
             break.x.by = 300,
             font.xtickslab = 15, 
             font.ytickslab = 15, 
             font.legend = 12,
             xlim = c(0, 1800),
             pval.coord = c(0, 0.05))
ggsurvc$table <- ggsurvc$table +
  ylab(NULL) + 
  xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsurvc
dev.off()


