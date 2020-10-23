
# Description -------------------------------------------------------------

# Figure 2: The B8T signature stratifies patients with TMB (high) and TMB (low)
# tumors into cohorts with different survival


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
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")))

# Bin TMB -----------------------------------------------------------------

breaks <- quantile(col_data$FMOne.mutation.burden.per.MB, probs = seq(0, 1, 1/2), na.rm = T)

col_data$tmb_bins <- cut(col_data$FMOne.mutation.burden.per.MB, breaks, labels = c("Lo", "Hi"))

col_data <- col_data %>% 
        unite(tmb_b8t, tmb_bins, b_bin, t_bin, remove = F) %>% 
        mutate(tmb_b8t = factor(tmb_b8t, levels = c("Hi_hi_hi", "Hi_lo_hi", "Hi_lo_lo", "Hi_hi_lo", "Lo_hi_hi", "Lo_lo_hi", "Lo_lo_lo", "Lo_hi_lo")))


# Fig 2a: TMB Bins vs Survival --------------------------------------------

fig_2a <- col_data %>% 
        filter(!is.na(tmb_bins)) 

png(filename = "./figures/fig_2/fig_2a.png", width = 5, height = 3, units = "in", res = 288)
survfit(Surv(os, censOS) ~ tmb_bins, data = fig_2a) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95), 
                   xlab = "Overall Survival (Months)",
                   legend.labs = c("Lo", "Hi"), 
                   legend.title = "TMB",
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 2b: TMB Bins across B8T ---------------------------------------------

fig_2b <- col_data %>% 
        filter(!is.na(tmb_bins)) 

labels <- c("TMB Lo", "TMB Hi")
names(labels) <- c("Lo", "Hi")

ggplot(fig_2b, aes(x = b_cell, y = cd8_t_eff, color = log2(FMOne.mutation.burden.per.MB + 1))) + 
        geom_point(size = 7, alpha = 0.7) + 
        scale_color_viridis_c(end = 0.95) + 
        facet_wrap(~tmb_bins, labeller = labeller(tmb_bins = labels)) + 
        theme_minimal() + 
        labs(x = "B-cell Signature", y = "CD8+ T-cell Signature", color = "Log2(TMB-per-MB + 1)") +
        geom_vline(xintercept = 0, size = 2, alpha = 0.5) +
        geom_hline(yintercept = 0, size = 2, alpha = 0.5) + 
        theme(text = element_text(size = 15),
              legend.position = "top",
              legend.key.size = unit(.3, "in"),
              legend.text = element_text(size = 15),
              legend.title = element_text(vjust = .85),
              panel.grid = element_blank(),
              panel.spacing = unit(2, "lines"),
              panel.background = element_rect(fill = "gray95", color = NA), 
              strip.text = element_text(size = 15))
ggsave("./figures/fig_2/fig_2b_i.png", width = 5.5, height = 4)

table_fig_2b <- fig_2b %>%
        group_by(tmb_bins, t_bin, b_bin) %>% 
        summarize(num = n()) %>%
        group_by(tmb_bins) %>% 
        mutate(prop = num/sum(num))

ggplot(table_fig_2b, aes(x = b_bin, y = t_bin, color = as.factor(rep(1:2, each = 4)))) + 
        geom_point(aes(size = num)) + 
        scale_color_viridis_d(begin = 0.1, end = 0.5) + 
        facet_wrap(~tmb_bins, labeller = labeller(tmb_bins = labels)) + 
        theme_minimal() + 
        scale_size_area(max_size = 40) + 
        geom_text(aes(label = round(prop, 2)), color = "white",  size = 7) +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank(),
              panel.spacing = unit(2, "lines"),
              panel.background = element_rect(fill = "gray95", color = NA), 
              strip.text = element_text(size = 15)) +
        labs(x = "B-cell Signature", y = "CD8+ T-cell Signature") +
        scale_x_discrete(labels = c("Lo", "Hi")) +
        scale_y_discrete(labels = c("Lo", "Hi"))
ggsave("./figures/fig_2/fig_2b_ii.png", width = 6, height = 3.5)


# Fig 2c: TMB Hi vs B8T Survival ------------------------------------------

fig_2c <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(tmb_bins == "Hi")

png(filename = "./figures/fig_2/fig_2c.png", width = 5, height = 3, units = "in", res = 288)
survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2c) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "B8T\n(TMB = Hi)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 2d: TMB Lo vs B8T Survival ------------------------------------------

fig_2d <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(tmb_bins == "Lo")

png(filename = "./figures/fig_2/fig_2d.png", width = 5, height = 3, units = "in", res = 288)
survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2d) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "B8T\n(TMB = Lo)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 2e: B8T Hi/Hi vs TMB Survival ---------------------------------------

fig_2e <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(b8t == "hi_hi")

png(filename = "./figures/fig_2/fig_2e.png", width = 5, height = 3, units = "in", res = 288)
survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2e) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1),
                   legend.title = "TMB\n(B8T = Hi/Hi)",
                   legend.labs = c("Hi", "Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 2f: B8T Lo/Lo vs TMB Survival ---------------------------------------

fig_2f <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(b8t == "lo_lo")

png(filename = "./figures/fig_2/fig_2f.png", width = 5, height = 3, units = "in", res = 288)
survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2f) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1),
                   legend.title = "TMB\n(B8T = Lo/Lo)",
                   legend.labs = c("Hi", "Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 2g: Opposing effects of TMB and B8T ---------------------------------

fig_2g <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(tmb_b8t %in% c("Lo_hi_hi", "Hi_lo_lo"))

png(filename = "./figures/fig_2/fig_2g.png", width = 5, height = 3, units = "in", res = 288)
survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2g) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(2, end = 0.95, direction = -1),
                   legend.title = "TMB/B8T",
                   legend.labs = c("Hi/Lo/Lo", "Lo/Hi/Hi"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()
