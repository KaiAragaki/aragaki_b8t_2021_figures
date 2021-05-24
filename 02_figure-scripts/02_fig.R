
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
                   title = "A.",
                   xlab = "Overall Survival (Months)",
                   legend.labs = c("Hi", "Lo"), 
                   legend.title = "TMB",
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
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
        labs(x = "B-cell Signature", y = "CD8+ T-cell Signature", color = "Log2(TMB-per-MB + 1)", title = "B.") +
        geom_vline(xintercept = 0, size = 2, alpha = 0.5) +
        geom_hline(yintercept = 0, size = 2, alpha = 0.5) + 
        coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
        theme(text = element_text(size = 15),
              legend.position = c(.45, 1.43),
              legend.direction = "horizontal",
              legend.key.size = unit(.3, "in"),
              legend.text = element_text(size = 15),
              legend.title = element_text(vjust = .85),
              panel.grid = element_blank(),
              plot.title.position = "plot",
              plot.title = element_text(face = "bold", vjust = 12),
              plot.margin = unit(c(17, 5, 1, 1), "mm"),
              panel.spacing = unit(2, "lines"),
              panel.background = element_rect(fill = "gray95", color = NA), 
              strip.text = element_text(size = 15))
ggsave("./figures/fig_2/fig_2b.png", width = 5.5, height = 4)


# Fig 2c: TMB Hi vs B8T Survival ------------------------------------------

fig_2c <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(tmb_bins == "Hi")

png(filename = "./figures/fig_2/fig_2c.png", width = 5.5, height = 6, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2c) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "B8T\n(TMB = Hi)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   title = "C.",
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   pval.coord = c(0, 0.05),
                   legend = "right")
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
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

png(filename = "./figures/fig_2/fig_2d.png", width = 5.5, height = 6, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2d) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "B8T\n(TMB = Lo)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"),
                   xlab = "Overall Survival (Months)", 
                   title = "D.",
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
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
                   title = "E.",
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
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
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig 2f: Opposing effects of TMB and B8T ---------------------------------

fig_2f <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        filter(tmb_b8t %in% c("Lo_hi_hi", "Hi_lo_lo"))

png(filename = "./figures/fig_2/fig_2f.png", width = 5.7, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_2f) %>% 
        ggsurvplot(pval = round(mt["Lo_hi_hi", "Hi_lo_lo"], digits = 4),
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1),
                   legend.title = "TMB/B8T",
                   legend.labs = c("Hi/Lo/Lo", "Lo/Hi/Hi"),
                   xlab = "Overall Survival (Months)", 
                   title = "F.",
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
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
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig 2g_i: Opposing effects of TMB and B8T with B8T/TMB Hi/Hi/Hi --------------

fig_2g_i <- col_data %>% 
        dplyr::filter(!is.na(tmb_bins)) %>%
        dplyr::mutate(group = case_when(tmb_bins == "Hi" & b8t == "hi_hi" ~ "TMB Hi, B8T Hi/Hi",
                                        tmb_bins == "Hi" ~ "TMB Hi, Non-Hi/Hi B8T",
                                        tmb_bins == "Lo" & b8t == "hi_hi" ~ "TMB Lo, B8T Hi/Hi",
                                        TRUE ~ NA_character_))
        
png(filename = "./figures/fig_2/fig_2g.png", width = 6.2, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ group, data = fig_2g_i) %>% 
        ggsurvplot(risk.table = T, pval = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(3, end = 0.95, direction = -1),
                   legend.title = "TMB; B8T",
                   legend.labs = c("Hi; Hi/Hi", "Hi; Non-Hi/Hi", "Lo; Hi/Hi"),
                   xlab = "Overall Survival (Months)", 
                   title = "G.",
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   pval.coord = c(0, 0.05),
                   legend = "right")
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()
