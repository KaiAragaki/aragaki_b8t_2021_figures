
# Description -------------------------------------------------------------

# Figure 2: The B8T signature stratifies patients with TMB (high) and TMB (low)
# tumors into cohorts with different OS in response to anti-PD-L1 therapy.


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
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               tmb_bins = if_else(FMOne.mutation.burden.per.MB < 10, "Lo", "Hi"),
               tmb_bins = factor(tmb_bins, levels = c("Hi", "Lo"))) %>% 
        unite(tmb_b8t, tmb_bins, b_bin, t_bin, remove = F) %>% 
        mutate(tmb_b8t = factor(tmb_b8t, levels = c("Hi_hi_hi", "Hi_lo_hi", "Hi_lo_lo", "Hi_hi_lo", "Lo_hi_hi", "Lo_lo_hi", "Lo_lo_lo", "Lo_hi_lo")))


# Fig 2a: TMB Bins vs Survival --------------------------------------------

png(filename = "./figures/fig_2/fig_2a.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_bins, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1), 
                   xlab = "Overall Survival (Months)",
                   legend.labs = c("Hi", "Lo"), 
                   legend.title = "TMB",
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig 2b: TMB Bins across B8T ---------------------------------------------

fig_2b <- col_data %>% 
        filter(!is.na(tmb_bins)) 

labels <- c("TMB Lo", "TMB Hi")
names(labels) <- c("Lo", "Hi")

ggplot(fig_2b, aes(x = b_cell, y = cd8_rose, color = log2(FMOne.mutation.burden.per.MB + 1))) + 
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
ggsave("./figures/fig_2/fig_2b.png", width = 5.5, height = 4)


# Fig 2c: TMB Hi vs B8T Survival ------------------------------------------

fig_2c <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(tmb_bins == "Hi")

png(filename = "./figures/fig_2/fig_2c.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2c) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "B8T\n(TMB = Hi)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())

ggsurv
dev.off()


# Fig 2d: TMB Lo vs B8T Survival ------------------------------------------

fig_2d <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(tmb_bins == "Lo")

png(filename = "./figures/fig_2/fig_2d.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2d) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "B8T\n(TMB = Lo)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())

ggsurv
dev.off()


# Fig 2e: B8T Hi/Hi vs TMB Survival ---------------------------------------

fig_2e <- col_data %>% 
        filter(!is.na(tmb_bins))

mt <- pairwise_survdiff(Surv(os, censOS) ~ tmb_b8t, data = fig_2e)$p.value

fig_2e <- fig_2e %>%  
        filter(b8t == "hi_hi")

png(filename = "./figures/fig_2/fig_2e.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2e) %>% 
        ggsurvplot(pval = round(mt["Lo_hi_hi", "Hi_hi_hi"], digits = 4), 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1),
                   legend.title = "TMB\n(B8T = Hi/Hi)",
                   legend.labs = c("Hi", "Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())

ggsurv
dev.off()


# Fig 2f: Opposing effects of TMB and B8T ---------------------------------

fig_2f <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(tmb_b8t %in% c("Lo_hi_hi", "Hi_lo_lo"))

png(filename = "./figures/fig_2/fig_2f.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2f) %>% 
        ggsurvplot(pval = round(mt["Lo_hi_hi", "Hi_lo_lo"], digits = 4),
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1),
                   legend.title = "TMB/B8T",
                   legend.labs = c("Hi/Lo/Lo", "Lo/Hi/Hi"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())

ggsurv
dev.off()


# Fig 2h: B8T Lo/Lo vs TMB Survival ---------------------------------------


fig_2h <- col_data %>% 
        unite(sex_tmb, gender, tmb_bins) %>% 
        mutate(sex_tmb = factor(sex_tmb, levels = c("male_Hi", "male_Lo", "female_Hi", "female_Lo")))

png(filename = "./figures/fig_2/fig_2h.png", width = 5, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ sex_tmb, data = fig_2h) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "Sex/TMB",
                   legend.labs = c("M/Hi", "M/Lo", "F/Hi", "F/Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()
