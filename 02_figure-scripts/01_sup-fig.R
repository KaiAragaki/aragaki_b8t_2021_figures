
# Description -------------------------------------------------------------

# Supplemental Figure 1. Patients with complete responses have higher CD8TGS


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(gt)
library(broom)
library(viridis)
library(DESeq2)

starify <- function(value) {
        if (value > 0.05) {
                "NS"
        } else if (value > 0.01) {
                "*"
        } else if (value > 0.001) {
                "**"
        } else {
                "***"
        }
}

compare_bar <- function(x, xend, y) {
        annotate("segment", x = x, xend = xend, y = y, yend = y, size = 1.5)
}


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")


# Extract colData ---------------------------------------------------------

col_data <- colData(imvigor) %>% 
        as_tibble() %>% 
        mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo")) %>% 
        mutate(Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>%
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")))


# Fig S1B: BCGS vs Response -----------------------------------------------

fig_s1a <- col_data %>% 
        filter(Best.Confirmed.Overall.Response != "NE")

# Stats

fig_s1a_stat <- 
        fig_s1a %>% 
        pivot_wider(names_from = Best.Confirmed.Overall.Response, values_from = b_cell)

cr_pr <- 
        t.test(fig_s1a_stat$CR, fig_s1a_stat$PR) %>% 
        tidy()
cr_sd <-        
        t.test(fig_s1a_stat$CR, fig_s1a_stat$SD) %>% 
        tidy()
cr_pd <- 
        t.test(fig_s1a_stat$CR, fig_s1a_stat$PD) %>% 
        tidy()

table_fig_s1a <- 
        table(fig_s1a$Best.Confirmed.Overall.Response) %>%
        as.data.frame() %>% 
        column_to_rownames("Var1") %>% 
        t() %>% 
        as.data.frame()

ggplot(fig_s1a, aes(x = Best.Confirmed.Overall.Response, 
                   y = b_cell, 
                   color = Best.Confirmed.Overall.Response, 
                   fill = Best.Confirmed.Overall.Response)) +
        geom_boxplot(alpha = 1, color = "black", lwd = .75) + 
        geom_jitter(height = 0, size = 2, alpha = 0.7, width = 0.2, color = "black") +
        theme_minimal() +
        scale_fill_viridis_d(end = 0.95) + 
        annotation_custom(text_grob(label = starify(cr_pr$p.value), 
                                    hjust = 0, 
                                    size = 30), 
                          ymin = 1.15, ymax = 1.15, xmin = 1.43, xmax = 1.43) +
        annotation_custom(text_grob(label = starify(cr_sd$p.value), 
                                    hjust = 0, 
                                    size = 30), 
                          ymin = 1.5, ymax = 1.5, xmin = 1.83, xmax = 1.83) +
        annotation_custom(text_grob(label = starify(cr_sd$p.value), 
                                    hjust = 0, 
                                    size = 30), 
                          ymin = 1.8, ymax = 1.8, xmin = 2.33, xmax = 2.33) +
        compare_bar(1, 2, 1.1) +
        compare_bar(1, 3, 1.4) +
        compare_bar(1, 4, 1.7) +
        geom_text(aes(x = 1, y = 2.1, label = paste0("(", CR, ")")), data = table_fig_s1a, size = 6, inherit.aes = F) + 
        geom_text(aes(x = 2, y = 2.1, label = paste0("(", PR, ")")), data = table_fig_s1a, size = 6, inherit.aes = F) + 
        geom_text(aes(x = 3, y = 2.1, label = paste0("(", SD, ")")), data = table_fig_s1a, size = 6, inherit.aes = F) + 
        geom_text(aes(x = 4, y = 2.1, label = paste0("(", PD, ")")), data = table_fig_s1a, size = 6, inherit.aes = F) + 
        theme(text = element_text(size = 25),
              legend.position = "none",
              panel.grid = element_blank(),
              axis.line.x.bottom = element_line(color = "black", size = 1),
              axis.line.y.left = element_line(color = "black", size = 1)) + 
        labs(x = NULL, y = "B-cell Signature") +
        coord_cartesian(ylim = c(-1, 1), clip = "off") + 
        theme(plot.margin = unit(c(12, 1, 1, 1), "lines"))

ggsave("./figures/fig_s1/fig_s1a.png", width = 5, height = 7)


# Fig S1B: CD8TGS vs Response ---------------------------------------------

fig_s1b <- col_data %>% 
        filter(Best.Confirmed.Overall.Response != "NE")

# Stats

fig_s1b_stat <- 
        fig_s1b %>% 
        pivot_wider(names_from = Best.Confirmed.Overall.Response, 
                    values_from = cd8_rose)

cr_pr <- 
        t.test(fig_s1b_stat$CR, 
               fig_s1b_stat$PR) %>% 
        tidy()
cr_sd <-        
        t.test(fig_s1b_stat$CR, 
               fig_s1b_stat$SD) %>% 
        tidy()
cr_pd <- 
        t.test(fig_s1b_stat$CR, 
               fig_s1b_stat$PD) %>% 
        tidy()

table_fig_s1b <- 
        table(fig_s1b$Best.Confirmed.Overall.Response) %>%
        as.data.frame() %>% 
        column_to_rownames("Var1") %>% 
        t() %>% 
        as.data.frame()

ggplot(fig_s1b, aes(x = Best.Confirmed.Overall.Response, 
                   y = cd8_rose, 
                   color = Best.Confirmed.Overall.Response, 
                   fill = Best.Confirmed.Overall.Response)) +
        geom_boxplot(alpha = 1, color = "black", lwd = .75) + 
        geom_jitter(height = 0, size = 2, alpha = 0.7, width = 0.2, color = "black") +
        theme_minimal() +
        scale_fill_viridis_d(end = 0.95) + 
        annotation_custom(text_grob(label = starify(cr_pr$p.value), 
                                    hjust = 0, 
                                    size = 30), 
                          ymin = 1.15, ymax = 1.15, xmin = 1.33, xmax = 1.33) +
        annotation_custom(text_grob(label = starify(cr_sd$p.value), 
                                    hjust = 0, 
                                    size = 30), 
                          ymin = 1.5, ymax = 1.5, xmin = 1.83, xmax = 1.83) +
        annotation_custom(text_grob(label = starify(cr_sd$p.value), 
                                    hjust = 0, 
                                    size = 30), 
                          ymin = 1.8, ymax = 1.8, xmin = 2.3, xmax = 2.3) +
        compare_bar(1, 2, 1.1) +
        compare_bar(1, 3, 1.4) +
        compare_bar(1, 4, 1.7) +
        geom_text(aes(x = 1, y = 2.1, label = paste0("(", CR, ")")), data = table_fig_s1b, size = 6, inherit.aes = F) + 
        geom_text(aes(x = 2, y = 2.1, label = paste0("(", PR, ")")), data = table_fig_s1b, size = 6, inherit.aes = F) + 
        geom_text(aes(x = 3, y = 2.1, label = paste0("(", SD, ")")), data = table_fig_s1b, size = 6, inherit.aes = F) + 
        geom_text(aes(x = 4, y = 2.1, label = paste0("(", PD, ")")), data = table_fig_s1b, size = 6, inherit.aes = F) + 
        theme(text = element_text(size = 25),
              legend.position = "none",
              panel.grid = element_blank(),
              axis.line.x.bottom = element_line(color = "black", size = 1),
              axis.line.y.left = element_line(color = "black", size = 1)) + 
        labs(x = NULL, y = "CD8+ T-cell Signature") +
        coord_cartesian(ylim = c(-1, 1), clip = "off") + 
        theme(plot.margin = unit(c(12, 1, 1, 1), "lines"))

ggsave("./figures/fig_s1/fig_s1b.png", width = 5, height = 7)


# Fig S1C: B8T vs Response -----------------------------------------------

fig_s1c <- col_data %>% 
        filter(Best.Confirmed.Overall.Response != "NE")

ggplot(fig_s1c, aes(b8t, fill = Best.Confirmed.Overall.Response)) + 
        geom_bar(position = "fill") + 
        scale_fill_viridis_d(end = 0.95, direction = -1) + 
        theme_minimal() + 
        theme(legend.position = "top",
              text = element_text(size = 15),
              panel.grid = element_blank()) +
        scale_x_discrete(labels = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo")) +
        labs(fill = "", x = "B-cell/CD8+ T-cell Signature", y = "Proportion")

ggsave("./figures/fig_s1/fig_s1c.png", width = 4, height = 3)

