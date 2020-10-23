
# Description -------------------------------------------------------------

# Figure 1: The B8T signature stratifies patients with metastatic bladder cancer
# into cohorts with different OS post-CPI


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(viridis)
library(DESeq2)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")


# Extract colData ---------------------------------------------------------

col_data <- colData(imvigor) %>% 
        as_tibble() %>% 
        mutate(t_bin = if_else(cd8_t_eff > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo")) %>% 
        mutate(Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE")))


# Fig 1a: Signature Distribution ------------------------------------------

ggplot(col_data, aes(b_cell, fill = 1)) +
        geom_density(alpha = 0.8, color = NA) + 
        scale_fill_viridis_c(end = .95) + 
        xlab("B-cell Signature") +
        ylab("Density") + 
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank())
ggsave("./figures/fig_1/fig_1a_i.png", width = 4, height = 3)

ggplot(col_data, aes(cd8_t_eff, fill = 1)) +
        geom_density(alpha = 0.8, color = NA) + 
        scale_fill_viridis_c(end = .95) + 
        xlab("CD8+ T-cell Signature") +
        ylab("Density") + 
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank())
ggsave("./figures/fig_1/fig_1a_ii.png", width = 4, height = 3)


# Fig 1b: Survival between Signature Modes --------------------------------

png(filename = "./figures/fig_1/fig_1b_i.png", width = 4, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ b_bin, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "B-cell Signature", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   pval.coord = c(0, 0.05))
dev.off()

png(filename = "./figures/fig_1/fig_1b_ii.png", width = 4, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ t_bin, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "CD8+ T-cell Signature", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 1c: Response vs Signatures Scatter ----------------------------------

fig_1c <- filter(col_data, Best.Confirmed.Overall.Response != "NE")

ggplot(fig_1c, aes(x = b_cell, y = cd8_t_eff, color = Best.Confirmed.Overall.Response)) + 
        scale_color_viridis_d(end = 0.90, direction = -1) + 
        geom_point(size = 4, alpha = 0.7) + 
        theme_minimal() + 
        labs(color = "Best Response", x = "B-cell Signature", y = "CD8+ T-cell Signature") +
        theme(text = element_text(size = 15),              
              panel.grid = element_blank())
ggsave("./figures/fig_1/fig_1c.png", width = 6, height = 4)


# Fig 1d: Survival vs B8T --------------------------------------------------

fig_1d <- col_data %>% 
        unite(b_t, b_bin, t_bin) %>%
        mutate(b_t = factor(b_t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")))

png(filename = "./figures/fig_1/fig_1d.png", width = 5.5, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ b_t, data = fig_1d) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"), 
                   legend.title = "B-cell/CD8+ Sig.", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12, 
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# B-Cell & T-Cell Vs Response Group Fig 1h -------------------------------

fig_1e <- col_data %>% 
        unite(b_t, b_bin, t_bin) %>%
        mutate(b_t = factor(b_t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo"))) %>% 
        filter(Best.Confirmed.Overall.Response != "NE")

ggplot(fig_1e, aes(b_t, fill = Best.Confirmed.Overall.Response)) + 
        geom_bar(position = "fill") + 
        scale_fill_viridis_d(end = 0.95, direction = -1) + 
        theme_minimal() + 
        theme(legend.position = "top",
              text = element_text(size = 15),
              panel.grid = element_blank()) +
        scale_x_discrete(labels = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo")) +
        labs(fill = "", x = "B-cell/CD8+ T-cell Signature", y = "Proportion")

ggsave("./figures/fig_1/fig_1e.png", width = 4, height = 3)
