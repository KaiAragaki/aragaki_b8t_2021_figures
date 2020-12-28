
# Description -------------------------------------------------------------

# Figure 1: The B8T signature stratifies patients with advanced urothelial 
# carcinoma into cohorts with different OS in response to anti-PD-L1 
# therapy


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
        mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
               b_bin = if_else(b_cell   > 0, "hi", "lo")) %>% 
        mutate(Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")))


# Fig 1a: Signature Distribution ------------------------------------------

ggplot(col_data, aes(b_cell, fill = 1)) +
        geom_density(alpha = 0.8, color = NA) + 
        scale_fill_viridis_c(end = .95) + 
        xlab("B-cell Signature") +
        ylab("Density") + 
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank()) +
        coord_cartesian(xlim = c(-1, 1))
ggsave("./figures/fig_1/fig_1a_i.png", width = 4, height = 3)

ggplot(col_data, aes(cd8_rose, fill = 1)) +
        geom_density(alpha = 0.8, color = NA) + 
        scale_fill_viridis_c(end = .95) + 
        xlab("CD8+ T-cell Signature") +
        ylab("Density") + 
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank()) +
        coord_cartesian(xlim = c(-1, 1))
ggsave("./figures/fig_1/fig_1a_ii.png", width = 4, height = 3)


# Fig 1b: Survival between Signature Modes --------------------------------

png(filename = "./figures/fig_1/fig_1b_i.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b_bin, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "B-cell Signature", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()

png(filename = "./figures/fig_1/fig_1b_ii.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ t_bin, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   legend.labs = c("High", "Low"), 
                   legend.title = "CD8+ T-cell Signature", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig 1c: Response vs Signatures Scatter ----------------------------------

fig_1c <- filter(col_data, Best.Confirmed.Overall.Response != "NE")

ggplot(fig_1c, aes(x = b_cell, y = cd8_rose, color = Best.Confirmed.Overall.Response)) + 
        scale_color_viridis_d(end = 0.90, direction = -1) + 
        geom_point(size = 2, alpha = 0.7) + 
        theme_minimal() + 
        labs(color = "Best Response", x = "B-cell Gene Signature", y = "CD8+ T-cell Gene Signature") +
        theme(text = element_text(size = 15),              
              panel.grid = element_blank())
ggsave("./figures/fig_1/fig_1c.png", width = 6, height = 4)


# Fig 1d: Survival vs B8T --------------------------------------------------

png(filename = "./figures/fig_1/fig_1d.png", width = 5.5, height = 6, units = "in", res = 288)
ggsurv <- 
        survfit(Surv(os, censOS) ~ b8t, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"), 
                   legend.title = "B8T", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12, 
                   legend = "top",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())

ggsurv
dev.off()


# Fig 1e: Survival vs B8T (No Platinum) ------------------------------------

no_plat <- filter(col_data, Received.platinum == "N")

png(filename = "./figures/fig_1/fig_1e.png", width = 5.5, height = 6, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b8t, data = no_plat) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   legend.title = "B8T\n(No Platinum)", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())

ggsurv
dev.off()

# Fig 1f: Survival vs B8T (Platinum) --------------------------------------

plat <- filter(col_data, Received.platinum == "Y")

png(filename = "./figures/fig_1/fig_1f.png", width = 5.5, height = 6, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b8t, data = plat) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   legend.title = "B8T\n(Platinum-Treated)", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 12,
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())

ggsurv
dev.off()
