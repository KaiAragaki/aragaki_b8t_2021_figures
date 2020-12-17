
# Description -------------------------------------------------------------

# Supplemental Figure 5


# Prepare Workspace --------------------------------------------------

library(tidyverse)
library(DESeq2)
library(survival)
library(survminer)
library(viridis)
library(broom)
library(gridExtra)

# Read in Data ------------------------------------------------------------

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
View(col_data)
# Fig S5a: B-Cell or T-cell Signature Distribution ------------------------

fig_s5a_1=ggplot(col_data, aes(b_cell, fill = 1)) +
  geom_density(alpha = 0.8, color = NA) + 
  scale_fill_viridis_c(end = .95) + 
  xlab("B-cell Signature") +
  ylab("Density") + 
  labs(subtitle="A")+
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = "none",
        panel.grid = element_blank())

ggsave(fig_s5a_1,file="H:/Rusty_B8T/figures/fig_5/fig_s5a_i.png", width = 4, height = 3)


fig_s5a_2=ggplot(col_data, aes(cd8_t_eff, fill = 1)) +
  geom_density(alpha = 0.8, color = NA) + 
  scale_fill_viridis_c(end = .95) + 
  xlab("CD8+ T-cell Signature") +
  ylab("Density") + 
  
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = "none",
        panel.grid = element_blank())

ggsave(fig_s5a_2,file="H:/Rusty_B8T/figures/fig_5/fig_s5a_ii.png", width = 4, height = 3)

ggsave(grid.arrange(fig_s5a_1,fig_s5a_2, nrow = 1), file="H:/Rusty_B8T/figures/fig_5/Fig_s5a.png",, width=8, height=3)
# Fig S5b: B-Cell Signature vs Survival -----------------------------------

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_s5b.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b_bin, data = col_data) %>% 
  ggsurvplot(pval = T, 
             risk.table = T,
             risk.table.height = 0.22,
             risk.table.title = "No. at risk",
             palette = viridis(2, end = 0.95, direction = -1), 
             legend.labs = c("High", "Low"), 
             legend.title = "B-cell Signature", 
             xlab = "Overall Survival (Months)",
             legend = "right",
             subtitle="B",
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


# Fig S5c: T-Cell Signature vs Survival -----------------------------------

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_s5c.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ t_bin, data = col_data) %>% 
  ggsurvplot(pval = T, 
             risk.table = T,
             risk.table.height = 0.22,
             risk.table.title = "No. at risk",
             palette = viridis(2, end = 0.95, direction = -1), 
             legend.labs = c("High", "Low"), 
             legend.title = "T-cell Signature", 
             xlab = "Overall Survival (Months)",
             legend = "right",
             subtitle="C",
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

# Fig S5d: B-Cell or T-cell Signature Distribution across Sex -------------

labels <- c("Male", "Female")
names(labels) <- c("male", "female")

fig_s5d1=ggplot(col_data, aes(b_cell, fill = patient.gender, color = patient.gender)) +
  geom_density(color = NA) +
  scale_fill_discrete(type = c("#1D557D", "#EB9BC7")) +
  xlab("B-cell Signature") +
  ylab("Density") + 
  labs(fill = "Sex",subtitle="D") +
  facet_wrap(~patient.gender, labeller = labeller(patient.gender = labels)) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.text = element_text(size = 15),
        panel.spacing = unit(2, "lines"),
        panel.background = element_rect(fill = "gray95", color = NA))

ggsave(fig_s5d1,file="H:/Rusty_B8T/figures/fig_5/fig_s5d_i.png", width = 8, height = 3)

fig_s5d2=ggplot(col_data, aes(cd8_t_eff, fill = patient.gender, color = patient.gender)) +
  geom_density(color = NA) +
  scale_fill_discrete(type = c("#1D557D", "#EB9BC7")) +
  xlab("CD8+ T-cell Signature") +
  ylab("Density") + 
  labs(fill = "Sex") +
  
  facet_wrap(~patient.gender, labeller = labeller(patient.gender = labels)) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.text = element_text(size = 15),
        panel.spacing = unit(2, "lines"),
        panel.background = element_rect(fill = "gray95", color = NA))
  

ggsave(fig_s5d2,file="H:/Rusty_B8T/figures/fig_5/fig_s5d_ii.png", width = 8, height = 3)
ggsave(grid.arrange(fig_s5d1,fig_s5d2, nrow = 2),,file="H:/Rusty_B8T/figures/fig_5/fig_s5d.png", width = 8, height = 6)


# Adjust P-values for Multiple Testing ------------------------------------

col_data_sex_b8t <- col_data %>% 
  unite(b8t_s, b8t, patient.gender, sep = ".")
table(col_data_sex_b8t$b8t_s)
mt <- pairwise_survdiff(Surv(new_death, death_event) ~ b8t_s, data = col_data_sex_b8t, p.adjust.method = "none") %>% 
  tidy() %>% 
  separate(group1, c("b8t_1", "sex_1"), sep = "\\.") %>% 
  separate(group2, c("b8t_2", "sex_2"), sep = "\\.") %>% 
  filter(b8t_1 == b8t_2) %>% 
  mutate(adj = p.adjust(p.value, "fdr"))
View(mt)
str(mt)
?pairwise_survdiff

# Compare Hi-Lo vs All others ----------------------------------------

col_data_2 <- col_data %>% 
  mutate(is_hi_lo = if_else(b8t == "hi_lo", T, F))

col_data_2_male <- filter(col_data_2, patient.gender == "male")

ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender + is_hi_lo, data = col_data_2_male) %>% 
  ggsurvplot(pval = T, 
             risk.table = T, surv.median.line = "hv",
             risk.table.height = 0.22,
             risk.table.title = "No. at risk",
             legend = "right",
             xlab = "Overall Survival (Months)", 
             xscale = 30,
             break.x.by = 300,
             font.xtickslab = 15, 
             font.ytickslab = 15, 
             font.legend = 10,
             xlim = c(0, 1800),
             pval.coord = c(0, 0.05))

col_data_2_female <- filter(col_data_2, patient.gender == "female")

ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender + is_hi_lo, data = col_data_2_female) %>% 
  ggsurvplot(pval = T, 
             risk.table = T, surv.median.line = "hv",
             risk.table.height = 0.22,
             risk.table.title = "No. at risk",
             legend = "right",
             xlab = "Overall Survival (Months)", 
             xscale = 30,
             break.x.by = 300,
             font.xtickslab = 15, 
             font.ytickslab = 15, 
             font.legend = 10,
             xlim = c(0, 1800),
             pval.coord = c(0, 0.05))



# Fig S5e Compare Sex Outcomes Between B8T Hi/Hi --------------------------

fig_s5e <- col_data %>% 
  filter(b8t == "hi_hi")

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_s5e.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender, data = fig_s5e) %>% 
  ggsurvplot(pval = T,#round(mt$adj[mt$b8t_1 == "hi_hi"], 2), #
             risk.table = T,
             risk.table.height = 0.22,
             risk.table.title = "No. at risk",
             palette = c("#1D557D", "#EB9BC7"),
             legend.labs = c("Male", "Female"), 
             legend.title = "Sex\n(B8T Hi/Hi)", 
             legend = "right",
             xlab = "Overall Survival (Months)",
             subtitle = "E",
             xscale = 30,
             break.x.by = 300,
             font.xtickslab = 15, 
             font.ytickslab = 15, 
             font.legend = 10,
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


# Fig S5f Compare Sex Outcomes Between B8T Hi/Lo --------------------------

fig_s5f <- col_data %>% 
  filter(b8t == "hi_lo")

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_s5f.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender, data = fig_s5f) %>% 
  ggsurvplot(pval =T,# round(mt$adj[mt$b8t_1 == "hi_lo"], 2), #
             risk.table = T,
             risk.table.height = 0.22,
             risk.table.title = "No. at risk",
             palette = c("#1D557D", "#EB9BC7"),
             legend.labs = c("Male", "Female"), 
             legend.title = "Sex\n(B8T Hi/Lo)", 
             legend = "right",
             xlab = "Overall Survival (Months)", 
             subtitle = "F",
             xscale = 30,
             break.x.by = 300,
             font.xtickslab = 15, 
             font.ytickslab = 15, 
             font.legend = 10,
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


# Fig S5g Compare Sex Outcomes Between B8T Lo/Hi --------------------------

fig_s5g <- col_data %>% 
  filter(b8t == "lo_hi")

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_s5g.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender, data = fig_s5g) %>% 
  ggsurvplot(pval =T,# round(mt$adj[mt$b8t_1 == "lo_hi"], 2),# 
             risk.table = T,
             risk.table.height = 0.22,
             risk.table.title = "No. at risk",
             palette = c("#1D557D", "#EB9BC7"),
             legend.labs = c("Male", "Female"), 
             legend.title = "Sex\n(B8T Lo/Hi)", 
             legend = "right",
             xlab = "Overall Survival (Months)", 
             subtitle = "G",
             xscale = 30,
             break.x.by = 300,
             font.xtickslab = 15,#8,# 
             font.ytickslab = 15, 
             font.legend = 10,
             #font.risk.table=3,#
             xlim = c(0, 5000),
             pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
  ylab(NULL) + 
  xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsurv
dev.off()


# Fig S5h Compare Sex Outcomes Between B8T Lo/Lo --------------------------

fig_s5h <- col_data %>% 
  filter(b8t == "lo_lo")

png(filename = "H:/Rusty_B8T/figures/fig_5/fig_s5h.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender, data = fig_s5h) %>% 
  ggsurvplot(pval =T,# round(mt$adj[mt$b8t_1 == "lo_lo"], 2), #
             risk.table = T,
             risk.table.height = 0.22,
             risk.table.title = "No. at risk",
             palette = c("#1D557D", "#EB9BC7"),
             legend.labs = c("Male", "Female"), 
             legend.title = "Sex\n(B8T Lo/Lo)", 
             legend = "right",
             xlab = "Overall Survival (Months)", 
             subtitle = "H",
             xscale = 30,
             break.x.by = 300,
             font.xtickslab = 15, 
             font.ytickslab = 15, 
             font.legend = 10,
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


