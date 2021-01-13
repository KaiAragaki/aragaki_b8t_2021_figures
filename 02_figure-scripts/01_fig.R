
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

col_data_long <- 
        pivot_longer(col_data, cols = c(b_cell, cd8_rose)) %>% 
        mutate(name = factor(name, levels = c("b_cell", "cd8_rose"), 
                             labels = c("B-cell Signature", "CD8+ T-cell Signature")))

ggplot(col_data_long, aes(value, fill = 1)) +
        facet_grid(~name) + 
        geom_density(alpha = 0.8, color = NA) + 
        scale_fill_viridis_c(end = .95) + 
        labs(x = "B-cell Signature", y = "Density", title = "A.") +
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank(), 
              plot.title.position = "plot",
              plot.title = element_text(face = "bold"),
              axis.title.x = element_blank(),
              panel.spacing = unit(4, "lines"),
              plot.margin = unit(c(1, 10, 1, 1), "mm")) +
        coord_cartesian(xlim = c(-1, 1), ylim = c(0, NA), expand = F)
ggsave("./figures/fig_1/fig_1a.png", width = 7, height = 3)


# Fig 1b: Survival between Signature Modes --------------------------------

png(filename = "./figures/fig_1/fig_1b_i.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b_bin, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1),
                   title = "B.",
                   legend.labs = c("High", "Low"), 
                   legend.title = "B-cell Signature", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 20, 
                   font.ytickslab = 20, 
                   font.legend = 20,
                   fontsize = 7,
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27))
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
                   title = " ",
                   legend.labs = c("High", "Low"), 
                   legend.title = "CD8+ T-cell Sig.", 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 20, 
                   font.ytickslab = 20, 
                   font.legend = 20,
                   fontsize = 7,
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27))
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
        labs(color = "Best Response", x = "B-cell Signature", 
             y = "CD8+ T-cell Signature", title = "C.") +
        coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = F) + 
        theme(text = element_text(size = 12),              
              panel.grid = element_blank(),
              plot.title.position = "plot",
              plot.title = element_text(face = "bold"),
              legend.position = "top",
              plot.margin = unit(c(1, 15, 1, 3), "mm"))
ggsave("./figures/fig_1/fig_1c.png", width = 4.3, height = 4)


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
                   title = "D.",
                   font.xtickslab = 17, 
                   font.ytickslab = 17, 
                   font.legend = 17,
                   fontsize = 6,
                   legend = "top",
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              plot.margin = unit(c(1, 10, 1, 1), "mm"))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.margin = unit(c(1, 10, 1, 1), "mm"))

ggsurv
dev.off()
