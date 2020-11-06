
# Description -------------------------------------------------------------

# Supplemental Figure 5


# Prepare Workspace --------------------------------------------------

library(tidyverse)
library(DESeq2)
library(survival)
library(survminer)
library(viridis)
library(broom)


# Read in Data ------------------------------------------------------------

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

# Fig S5a: B-Cell or T-cell Signature Distribution ------------------------

ggplot(col_data, aes(b_cell, fill = 1)) +
        geom_density(alpha = 0.8, color = NA) + 
        scale_fill_viridis_c(end = .95) + 
        xlab("B-cell Signature") +
        ylab("Density") + 
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank())

ggsave("./figures/fig_s5/fig_s5a_i.png", width = 4, height = 3)


ggplot(col_data, aes(cd8_t_eff, fill = 1)) +
        geom_density(alpha = 0.8, color = NA) + 
        scale_fill_viridis_c(end = .95) + 
        xlab("CD8+ T-cell Signature") +
        ylab("Density") + 
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank())

ggsave("./figures/fig_s5/fig_s5a_ii.png", width = 4, height = 3)


# Fig S5b: B-Cell Signature vs Survival -----------------------------------

png(filename = "./figures/fig_s5/fig_s5b.png", width = 5.5, height = 5.5, units = "in", res = 288)
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

png(filename = "./figures/fig_s5/fig_s5c.png", width = 5.5, height = 5.5, units = "in", res = 288)
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

ggplot(col_data, aes(b_cell, fill = patient.gender, color = patient.gender)) +
        geom_density(color = NA) +
        scale_fill_discrete(type = c("#1D557D", "#EB9BC7")) +
        xlab("B-cell Signature") +
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

ggsave("./figures/fig_s5/fig_s5d_i.png", width = 6, height = 3)

ggplot(col_data, aes(cd8_t_eff, fill = patient.gender, color = patient.gender)) +
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

ggsave("./figures/fig_s5/fig_s5d_ii.png", width = 6, height = 3)



# Adjust P-values for Multiple Testing ------------------------------------

col_data_sex_b8t <- col_data %>% 
        unite(b8t_s, b8t, patient.gender, sep = ".")

mt <- pairwise_survdiff(Surv(new_death, death_event) ~ b8t_s, data = col_data_sex_b8t, p.adjust.method = "none") %>% 
        tidy() %>% 
        separate(group1, c("b8t_1", "sex_1"), sep = "\\.") %>% 
        separate(group2, c("b8t_2", "sex_2"), sep = "\\.") %>% 
        filter(b8t_1 == b8t_2) %>% 
        mutate(adj = p.adjust(p.value, "fdr"))


# Fig S5e Compare Sex Outcomes Between B8T Hi/Hi --------------------------

fig_s5e <- col_data %>% 
        filter(b8t == "hi_hi")

png(filename = "./figures/fig_s5/fig_s5e.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender, data = fig_s5e) %>% 
        ggsurvplot(pval = round(mt$adj[mt$b8t_1 == "hi_hi"], 2), 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = c("#1D557D", "#EB9BC7"),
                   legend.labs = c("Male", "Female"), 
                   legend.title = "Sex\n(B8T Hi/Hi)", 
                   legend = "right",
                   xlab = "Overall Survival (Months)", 
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

png(filename = "./figures/fig_s5/fig_s5f.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender, data = fig_s5f) %>% 
        ggsurvplot(pval = round(mt$adj[mt$b8t_1 == "hi_lo"], 2), 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = c("#1D557D", "#EB9BC7"),
                   legend.labs = c("Male", "Female"), 
                   legend.title = "Sex\n(B8T Hi/Lo)", 
                   legend = "right",
                   xlab = "Overall Survival (Months)", 
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

png(filename = "./figures/fig_s5/fig_s5g.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender, data = fig_s5g) %>% 
        ggsurvplot(pval = round(mt$adj[mt$b8t_1 == "lo_hi"], 2), 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = c("#1D557D", "#EB9BC7"),
                   legend.labs = c("Male", "Female"), 
                   legend.title = "Sex\n(B8T Lo/Hi)", 
                   legend = "right",
                   xlab = "Overall Survival (Months)", 
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


# Fig S5h Compare Sex Outcomes Between B8T Lo/Lo --------------------------

fig_s5h <- col_data %>% 
        filter(b8t == "lo_lo")

png(filename = "./figures/fig_s5/fig_s5h.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ patient.gender, data = fig_s5h) %>% 
        ggsurvplot(pval = round(mt$adj[mt$b8t_1 == "lo_lo"], 2), 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = c("#1D557D", "#EB9BC7"),
                   legend.labs = c("Male", "Female"), 
                   legend.title = "Sex\n(B8T Lo/Lo)", 
                   legend = "right",
                   xlab = "Overall Survival (Months)", 
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


# Read in Melanoma Data ---------------------------------------------------

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


# Filter for Primary Tumors -----------------------------------------------

pri <- filter(col_data, patient.tumor_tissue_site == "primary tumor")


# Fig S5i: Female Primary Melanoma vs Survival across B8T -----------------

fig_s5i <- filter(pri, patient.gender == "female")

png(filename = "./figures/fig_s5/fig_s5i.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_spec, data = fig_s5i) %>% 
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
                   break.x.by = 365,
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


# Fig S5j: Male Primary Melanoma vs Survival across B8T -------------------

fig_s5j <- filter(pri, patient.gender == "male")

png(filename = "./figures/fig_s5/fig_s5j.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_spec, data = fig_s5j) %>% 
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
                   break.x.by = 365,
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   xlim = c(0, 1500),
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Filter for Regional/Distant Metastases ----------------------------------

reg_met <- filter(col_data, patient.tumor_tissue_site %in% c("regional cutaneous or subcutaneous tissue (includes satellite and in-transit metastasis)", "distant metastasis"))


# Fig S5k: Female Regional/Distant Metastases vs Survival across B8T ------

fig_s5k <- filter(reg_met, patient.gender == "female")

png(filename = "./figures/fig_s5/fig_s5k.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_spec, data = fig_s5k) %>% 
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


# Fig S5l: Male Regional/Distant Metastases vs Survival across B8T --------

fig_s5l <- filter(reg_met, patient.gender == "male")

png(filename = "./figures/fig_s5/fig_s5l.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(new_death, death_event) ~ b8t_spec, data = fig_s5l) %>% 
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
